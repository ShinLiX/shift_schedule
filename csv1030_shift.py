#!/usr/bin/env python3
import csv
from ortools.sat.python import cp_model

# ==================== 1) Scheduling Logic ====================

def compute_total_weekly_availability(employees_availability):
    """Return {employee: total # of available weekly slots (0..83)}."""
    return { e: len(avail) for e, avail in employees_availability.items() }

def is_contiguous_ones(bits):
    """Check if bits is exactly one contiguous block of 1s, or all zero."""
    ones_positions = [i for i,b in enumerate(bits) if b==1]
    if not ones_positions:
        return True
    first = ones_positions[0]
    last  = ones_positions[-1]
    return (last - first + 1) == len(ones_positions)

def generate_day_patterns(availability_mask, min_block=6, max_block=8):
    """
    For a single day of 12 hours => sum=0 => off,
    if sum(availability_mask)<6 => fill smaller block if contiguous or off,
    else => one contiguous block in [6..8] or off.
    Return patterns of length=13 => 12 bits + 1 bit for 'workedBit'.
    """
    n = sum(availability_mask)
    valid_pats = []
    for mask in range(1<<12):
        pat = [(mask >> i)&1 for i in range(12)]
        # Must be subset
        conflict=False
        for i in range(12):
            if pat[i]==1 and availability_mask[i]==0:
                conflict=True
                break
        if conflict:
            continue
        total_ones = sum(pat)
        worked_bit = 1 if total_ones>0 else 0

        if n< min_block:
            # If employee's daily availability is <6 => they can fill that entire smaller block (if contiguous) or off
            if total_ones==0:
                valid_pats.append(pat+[0])
            elif total_ones==n:
                # fill entire smaller block if contiguous
                same=True
                for i in range(12):
                    if pat[i]!= availability_mask[i]:
                        same=False
                        break
                if same and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
        else:
            # else => 6..8 hours or off
            if total_ones==0:
                valid_pats.append(pat+[0])
            else:
                if (min_block<=total_ones<=max_block) and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
    return valid_pats

def solve_scheduling(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    coverage_needed=2,
    desired_shifts=None
):
    """
    EXACT coverage=2, single shift/day approach:
      - If desired_shifts[e] in [1,2], forcibly total_shifts[e] = that many
      - If desired_shifts[e]>2 => we do a soft approach, never exceed desired_shifts
      - override min hours if availability < global_min or desired_shifts[e]*8< global_min
    Returns: (schedule, solver, status, total_shifts_vars)
      - schedule => {emp: [slots]}
      - solver => the CpSolver
      - status => solver status
      - total_shifts_vars => dict {emp: the IntVar for total shifts}
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    employees = list(employees_availability.keys())
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees}

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    # x[e,s]=1 => e works slot s
    x={}
    for e in employees:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage=2 => sum(x[e,s])=2
    for s in range(total_slots):
        model.Add( sum(x[(e,s)] for e in employees)== coverage_needed )

    # if not in availability => x[e,s]=0
    for e in employees:
        for s in range(total_slots):
            if s not in employees_availability[e]:
                model.Add(x[(e,s)]==0)

    # override min
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees:
        ds= desired_shifts[e]
        if (tot_avail[e]< global_min_hours) or (ds*8< global_min_hours):
            eff_min= 0
        else:
            eff_min= global_min_hours
        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add( weekly_hours[e]>= eff_min)
        model.Add( weekly_hours[e]<= global_max_hours)

    # forbidden pairs => x[e1,s] + x[e2,s] <=1
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees and emp2 in employees:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)]<=1)

    # day-based pattern => each day => patterns
    day_pattern_index={}
    day_worked={}
    all_patterns={}
    for e in employees:
        all_patterns[e]={}
        for d in range(days):
            avmask=[]
            startSlot= d*12
            for i in range(12):
                slot_id= startSlot + i
                avmask.append(1 if slot_id in employees_availability[e] else 0)
            valid_pats= generate_day_patterns(avmask,6,8)
            patIdx= model.NewIntVar(0,len(valid_pats)-1, f"patIndex_{e}_{d}")
            day_pattern_index[(e,d)] = patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    # build table constraints
    for e in employees:
        for d in range(days):
            day_vars= [ x[(e, d*12 + i)] for i in range(12) ]
            day_vars_plus= day_vars + [ day_worked[(e,d)], day_pattern_index[(e,d)] ]
            table=[]
            valid_pats= all_patterns[e][d]
            for idx, pat in enumerate(valid_pats):
                row= pat[:12] + [pat[12], idx]
                table.append(row)
            # constraint
            model.AddAllowedAssignments(day_vars_plus, table)

    # total_shifts[e]
    total_shifts={}
    for e in employees:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"total_shifts_{e}")
        model.Add( sum(sumDays)== total_shifts[e])

    # objective => if ds in [1,2], forced; else never exceed ds => minimize difference
    diffs=[]
    for e in employees:
        ds= desired_shifts[e]
        if ds in [1,2]:
            model.Add( total_shifts[e]== ds )
        else:
            model.Add( total_shifts[e]<= ds )
            diffUp= model.NewIntVar(0,7,f"diffUp_{e}")
            diffDown= model.NewIntVar(0,7,f"diffDown_{e}")
            model.Add( total_shifts[e]- ds <= diffUp )
            model.Add( ds- total_shifts[e] <= diffDown )
            diffs.append(diffUp)
            diffs.append(diffDown)
    model.Minimize( sum(diffs) )

    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        schedule={ e:[] for e in employees}
        for e in employees:
            for s in range(total_slots):
                if solver.Value(x[(e,s)])==1:
                    schedule[e].append(s)
        # Return also the dict of total_shifts[e] => the IntVar
        return schedule, solver, status, total_shifts
    else:
        return None, None, status, None

# ==================== 2) CSV Export (include second table: Name, WeeklyHours, DesiredShift, ActualShift) ====================

def export_schedule_to_csv(
    schedule, employees, solver, total_shifts_vars,
    desired_shifts,
    filename="schedule.csv"
):
    """
    We produce columns: Time, Monday-1, Monday-2, Tuesday-1, Tuesday-2, ... Sunday-1, Sunday-2
    Then 12 time rows. Then a blank row, then a second table with headers:
      Name, WeeklyHours, DesiredShifts, ActualShifts
    one row per employee => each cell is separate.
    """
    day_names = ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]
    days=7
    slots_per_day=12

    # Build day_hour => list_of_emps
    day_hour_emps = { d: {h: [] for h in range(slots_per_day)} for d in range(days)}
    for e in employees:
        for slot_id in schedule[e]:
            d= slot_id// slots_per_day
            h= slot_id% slots_per_day
            day_hour_emps[d][h].append(e)

    # occupant continuity => day_col_occ[d][col][hour]
    day_col_occ= { d: {1: ["" for _ in range(slots_per_day)],
                        2: ["" for _ in range(slots_per_day)]}
                   for d in range(days)}

    for d in range(days):
        occupant1=None
        occupant2=None
        for h in range(slots_per_day):
            assigned= sorted(day_hour_emps[d][h])
            # keep occupant1 if occupant1 in assigned
            if occupant1 in assigned:
                assigned.remove(occupant1)
            else:
                occupant1=None
            if occupant2 in assigned:
                assigned.remove(occupant2)
            else:
                occupant2=None

            if occupant1 is None and assigned:
                occupant1= assigned.pop(0)
            if occupant2 is None and assigned:
                occupant2= assigned.pop(0)

            day_col_occ[d][1][h]= occupant1 if occupant1 else ""
            day_col_occ[d][2][h]= occupant2 if occupant2 else ""

    # define blocks => occupant stays => only print occupant in the first row
    day_col_blocks= { d: {1:[], 2:[]} for d in range(days)}
    for d in range(days):
        for col in [1,2]:
            occupant_list= day_col_occ[d][col]
            h=0
            while h< slots_per_day:
                nameHere= occupant_list[h]
                if nameHere=="":
                    day_col_blocks[d][col].append((h,h,""))
                    h+=1
                else:
                    startH= h
                    while (h+1<slots_per_day) and occupant_list[h+1]== nameHere:
                        h+=1
                    endH= h
                    day_col_blocks[d][col].append((startH,endH,nameHere))
                    h+=1

    occupantPrint= { d: {col: [""]*slots_per_day for col in [1,2]} for d in range(days)}
    for d in range(days):
        for col in [1,2]:
            for (startH,endH,name) in day_col_blocks[d][col]:
                occupantPrint[d][col][startH]= name

    # build CSV rows => 12 times => 11..12, 12..13, ... 22..23
    headers=["Time"]
    for d in range(days):
        day_name= day_names[d]
        headers.append(f"{day_name}-1")
        headers.append(f"{day_name}-2")

    rows=[]
    for hour in range(slots_per_day):
        time_str= f"{11+hour}:00 - {12+hour}:00"
        row_cells=[time_str]
        for d in range(days):
            row_cells.append( occupantPrint[d][1][hour] )
            row_cells.append( occupantPrint[d][2][hour] )
        rows.append(row_cells)

    # Then a blank row
    rows.append([])

    # Then a second table with a header row => "Name", "WeeklyHours", "DesiredShifts", "ActualShifts"
    # We'll store that in columns 0..3, ignoring the day columns
    rows.append(["Name","WeeklyHours","DesiredShifts","ActualShifts"])
    for e in employees:
        totalH= len(schedule[e])
        wanted= desired_shifts[e]
        # actual => solver.Value(total_shifts_vars[e])
        got= solver.Value(total_shifts_vars[e])
        rows.append([ e, str(totalH), str(wanted), str(got) ])

    # now write
    with open(filename,"w", newline="", encoding="utf-8") as f:
        writer= csv.writer(f)
        writer.writerow(headers)
        for r in rows:
            writer.writerow(r)

    print(f"Exported schedule to '{filename}'.")
    print("First table => time rows x day columns. Then a blank row, then Name/WeeklyHours/DesiredShifts/ActualShifts table.")


# ==================== 3) Putting it all together (interactive) ====================

def user_input_schedule():
    print("=== SHIFT SCHEDULER: EXACT=2, single shift. If desired=1 or2 => forced, else soft approach.")
    print("We use days=1..7 => Monday..Sunday, hours in [11..23].")

    global_min= int(input("Enter global min weekly hours (16): "))
    global_max= int(input("Enter global max weekly hours (40): "))
    nemp= int(input("How many employees? "))

    employees_availability={}
    desired_shifts={}

    for i in range(nemp):
        name= input(f"\nEmployee#{i+1} name: ").strip()
        print("Enter day-based availability => 'day start end' or 'all' or 'done', day in [1..7].")
        print(" e.g. '1 11 17' => Monday from 11..17. 'all' => 7*12=84 slots.")
        avail=set()
        while True:
            line=input("> ").strip().lower()
            if line=="done":
                break
            if line=="all":
                avail= set(range(7*12))
                break
            parts=line.split()
            if len(parts)==3:
                try:
                    day= int(parts[0])
                    sh= int(parts[1])
                    eh= int(parts[2])
                    if 1<= day<=7 and 11<= sh<eh<=23:
                        # convert to slot_id => day_idx= day-1
                        day_idx= day-1
                        for hour in range(sh,eh):
                            slot_id= day_idx*12 + (hour-11)
                            avail.add(slot_id)
                    else:
                        print("Day or hour out of range.")
                except:
                    print("Invalid. 'day start end' or 'all' or 'done'")
            else:
                print("type 'day start end' or 'all' or 'done'")
        employees_availability[name]= avail
        ds= int(input(f"Enter {name}'s desired #shifts => if 1 or2 => forced: "))
        desired_shifts[name]= ds

    # forbidden pairs
    forbidden=[]
    print("\nEnter forbidden pairs(emp1 emp2) or 'none'/'done' to skip.")
    while True:
        line=input("Forbidden pair> ").strip()
        if line in ("none","done"):
            break
        parts=line.split()
        if len(parts)==2:
            forbidden.append((parts[0], parts[1]))
        else:
            print("type 'emp1 emp2' or 'none'/'done'")

    schedule, solver, status, total_shifts_vars = solve_scheduling(
        employees_availability,
        global_min_hours=global_min,
        global_max_hours=global_max,
        forbidden_pairs= forbidden,
        coverage_needed=2,
        desired_shifts= desired_shifts
    )
    if schedule is None:
        print(f"No solution found. status={status}")
    else:
        print(f"Solution found! CP-SAT status={status}")
        # show each employee's weekly hours
        for e in employees_availability.keys():
            total_hrs= len(schedule[e])
            print(f"  {e} => {total_hrs} hours, desired={desired_shifts[e]}, got_shifts={solver.Value(total_shifts_vars[e])}")

        ans= input("Export to CSV? y/n: ").strip().lower()
        if ans=="y":
            fname= input("Enter csv filename (default='schedule.csv'): ").strip()
            if not fname:
                fname="schedule.csv"
            export_schedule_to_csv(
                schedule,
                list(employees_availability.keys()),
                solver,
                total_shifts_vars,
                desired_shifts,
                filename=fname
            )
        else:
            print("Ok, not exporting CSV.")


def main():
    user_input_schedule()

if __name__=="__main__":
    main()