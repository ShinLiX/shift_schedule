#!/usr/bin/env python3
import tkinter as tk
from tkinter import messagebox, filedialog
import csv
from ortools.sat.python import cp_model

# ========================== 1) Core Scheduling Logic - Unchanged for Main ==========================
def compute_total_weekly_availability(employees_availability):
    return { e: len(avail) for e, avail in employees_availability.items() }

def is_contiguous_ones(bits):
    idxs = [i for i,b in enumerate(bits) if b==1]
    if not idxs:
        return True
    first = idxs[0]
    last  = idxs[-1]
    return (last - first + 1) == len(idxs)

def generate_day_patterns(availability_mask, min_block=6, max_block=8):
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
            if total_ones==0:
                valid_pats.append(pat+[0])
            elif total_ones==n:
                same=True
                for i in range(12):
                    if pat[i]!= availability_mask[i]:
                        same=False
                        break
                if same and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
        else:
            if total_ones==0:
                valid_pats.append(pat+[0])
            else:
                if (min_block<=total_ones<=max_block) and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
    return valid_pats

# ========================== 2) Main Solve: coverage=2 EXACT ==========================
def solve_scheduling_main(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    coverage_needed=2,
    desired_shifts=None
):
    """
    EXACT coverage=2 approach (unchanged from your original).
    Returns schedule or None if infeasible.
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

    x={}
    for e in employees:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage=2 => sum(x[e,s])=2
    for s in range(total_slots):
        model.Add(sum(x[(e,s)] for e in employees)== coverage_needed)

    # availability
    for e in employees:
        for s in range(total_slots):
            if s not in employees_availability[e]:
                model.Add(x[(e,s)]==0)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees:
        ds= desired_shifts[e]
        eff_min= 0 if (tot_avail[e]< global_min_hours or ds*8< global_min_hours) else global_min_hours
        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add(weekly_hours[e]>= eff_min)
        model.Add(weekly_hours[e]<= global_max_hours)

    # forbidden pairs
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees and emp2 in employees:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)]<=1)

    # day-based pattern
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

    for e in employees:
        for d in range(days):
            day_vars= [ x[(e,d*12 + i)] for i in range(12)]
            plus= day_vars + [day_worked[(e,d)], day_pattern_index[(e,d)]]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12]+ [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # total_shifts
    total_shifts={}
    diffs=[]
    for e in employees:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"total_shifts_{e}")
        model.Add(sum(sumDays)== total_shifts[e])

    for e in employees:
        ds= desired_shifts[e]
        if ds in [1,2]:
            model.Add(total_shifts[e]== ds)
        else:
            model.Add(total_shifts[e]<= ds)
            diffUp= model.NewIntVar(0,7,f"diffUp_{e}")
            diffDown= model.NewIntVar(0,7,f"diffDown_{e}")
            model.Add(total_shifts[e]- ds <= diffUp)
            model.Add(ds - total_shifts[e] <= diffDown)
            diffs.append(diffUp)
            diffs.append(diffDown)

    model.Minimize(sum(diffs))
    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    schedule= {e:[] for e in employees}
    for e in employees:
        for s in range(total_slots):
            if solver.Value(x[(e,s)])==1:
                schedule[e].append(s)
    return schedule, solver, status, total_shifts


# ========================== 3) Fallback Solve: coverage ≤2, maximize coverage ==========================
def solve_scheduling_fallback(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    desired_shifts=None
):
    """
    If main solve fails, we do partial coverage:
     coverage <=2 each hour, and we maximize sum of x[e,s].
    We keep the same single-shift constraints, min/max hours, forbidden pairs, desired shifts, etc.
    This yields partial coverage. Some slots might remain uncovered => NA in CSV.
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    employees= list(employees_availability.keys())
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees}

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    x={}
    for e in employees:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage <=2 => sum(x[e,s]) <= 2
    for s in range(total_slots):
        model.Add( sum(x[(e,s)] for e in employees ) <= 2 )

    # availability
    for e in employees:
        for s in range(total_slots):
            if s not in employees_availability[e]:
                model.Add(x[(e,s)]==0)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees:
        ds= desired_shifts[e]
        eff_min= 0 if (tot_avail[e]< global_min_hours or ds*8< global_min_hours) else global_min_hours
        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add(weekly_hours[e]>= eff_min)
        model.Add(weekly_hours[e]<= global_max_hours)

    # forbidden
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees and emp2 in employees:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)]<=1)

    # day-based single shift
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
            patIdx= model.NewIntVar(0,len(valid_pats)-1,f"patIndex_{e}_{d}")
            day_pattern_index[(e,d)]= patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    for e in employees:
        for d in range(days):
            day_vars= [x[(e,d*12 + i)] for i in range(12)]
            plus= day_vars + [ day_worked[(e,d)], day_pattern_index[(e,d)] ]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12] + [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # desired shifts logic (like main, we do min difference)
    # but let's do a simpler approach => either the same approach or do min difference. 
    # We'll do the same approach for consistency:

    total_shifts={}
    diffs=[]
    for e in employees:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"fallback_shifts_{e}")
        model.Add(sum(sumDays)== total_shifts[e])
        ds= desired_shifts[e]
        if ds in [1,2]:
            model.Add(total_shifts[e]== ds)
        else:
            model.Add(total_shifts[e]<= ds)
            diffUp= model.NewIntVar(0,7,f"fallback_diffUp_{e}")
            diffDown= model.NewIntVar(0,7,f"fallback_diffDown_{e}")
            model.Add(total_shifts[e]- ds <= diffUp )
            model.Add(ds - total_shifts[e] <= diffDown )
            diffs.append(diffUp)
            diffs.append(diffDown)

    # But now we want to maximize coverage => sum x[e,s]
    # We can define coverage_sum= sum(x[e,s] for all e,s).
    coverage_sum= model.NewIntVar(0, 7*12*2, "coverage_sum")  # up to 168 if 7*12=84 slots *2 coverage
    model.Add( coverage_sum == sum(x[(e,s)] for e in employees for s in range(total_slots)) )

    # We'll define an objective => maximize coverage_sum - sum(diffs)
    # or we can keep the "min sum(diffs)" approach. In your original code, you did min sum(diffs).
    # We want partial coverage primarily, so let's do a 2-part objective => either Weighted approach or simpler approach => 
    # We'll do a Weighted approach => e.g. model.Maximize( coverage_sum * 1000 - sum(diffs) ) 
    # so we first try to get maximum coverage, then also try to reduce shift differences. 
    # This is a heuristic. Another approach is to min sum(diffs) as sub-objective. 
    # For demonstration, let's do coverage_sum main priority, then minimize sum(diffs).
    # => coverage_sum in integer form, sum(diffs) also. We'll define coverage_sum as well. 
    # Weighted approach => coverage_sum is to be maximized => so we'll do negative diffs => coverage_sum * 1000 - sum(diffs).
    # You can tune the weighting.

    model.Maximize( coverage_sum * 1000 - sum(diffs) )

    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    schedule= { e:[] for e in employees}
    for e in employees:
        for s in range(total_slots):
            if solver.Value(x[(e,s)])==1:
                schedule[e].append(s)
    return schedule, solver, status, total_shifts


# ========================== 4) Export (unchanged) ==========================
def export_schedule_to_csv(
    schedule, employees, solver, total_shifts_vars, desired_shifts, filename="schedule.csv"
):
    day_names = ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]
    days=7
    slots_per_day=12

    day_hour_emps = { d: {h: [] for h in range(slots_per_day)} for d in range(days)}
    # We'll also track coverage_count so we know if coverage=2,1,0 => fill occupant columns or "NA"
    coverage_count= { (d,h): 0 for d in range(days) for h in range(slots_per_day) }

    for e in employees:
        for slot_id in schedule[e]:
            d= slot_id// slots_per_day
            h= slot_id% slots_per_day
            day_hour_emps[d][h].append(e)

    # occupant continuity
    # but we must handle coverage <2 => occupant1=someone if any, occupant2=someone if coverage=2, else "NA"
    # We'll do the same occupant1 occupant2 approach, but if occupant2 doesn't exist => occupant2=NA
    # If occupant1 also doesn't exist => occupant1=NA

    day_col_occ= { d: {1: ["" for _ in range(slots_per_day)],
                        2: ["" for _ in range(slots_per_day)]}
                   for d in range(days)}

    for d in range(days):
        occupant1=None
        occupant2=None
        for h in range(slots_per_day):
            assigned= sorted(day_hour_emps[d][h])
            coverage_count[(d,h)] = len(assigned)  # 0..2
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

    # Build main table
    headers=["Time"]
    for d in range(days):
        headers.append(f"{day_names[d]}-1")
        headers.append(f"{day_names[d]}-2")

    rows=[]
    for hour in range(slots_per_day):
        time_str= f"{11+hour}:00 - {12+hour}:00"
        row_cells=[time_str]
        for d in range(days):
            c1= occupantPrint[d][1][hour]
            c2= occupantPrint[d][2][hour]
            # if coverage_count[(d,h)] <2 => occupant2= NA
            # if coverage_count[(d,h)]==0 => occupant1= NA
            # We will override occupant if coverage is partial
            cov= coverage_count[(d,h)]
            if cov==0:
                c1="NA"
                c2="NA"
            elif cov==1:
                # occupant2 => "NA"
                if not c1 and c2:
                    # occupant1 is blank but occupant2 is filled => swap
                    c1, c2 = c2, "NA"
                elif not c2:
                    c2="NA"
            row_cells.append(c1)
            row_cells.append(c2)
        rows.append(row_cells)

    rows.append([])
    rows.append(["Name","WeeklyHours","DesiredShifts","ActualShifts"])
    for e in employees:
        totalH= len(schedule[e])
        wanted= desired_shifts[e]
        got= solver.Value(total_shifts_vars[e])
        rows.append([ e, str(totalH), str(wanted), str(got) ])

    with open(filename,"w", newline="", encoding="utf-8") as f:
        writer= csv.writer(f)
        writer.writerow(headers)
        for r in rows:
            writer.writerow(r)

    print(f"Exported schedule to '{filename}'.")


# ========================== 6) Single-Page Tkinter GUI with 2-Phase Solve ==========================
class SinglePageSchedulerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Employee Scheduler - Single Page")

        fm_global= tk.LabelFrame(root, text="Global Constraints", padx=5, pady=5)
        fm_global.pack(fill='x', padx=10, pady=5)

        tk.Label(fm_global, text="Global min hours:").grid(row=0, column=0, sticky='e')
        self.global_min_var= tk.StringVar(value="16")
        tk.Entry(fm_global, textvariable=self.global_min_var, width=5).grid(row=0, column=1, padx=5)

        tk.Label(fm_global, text="Global max hours:").grid(row=0, column=2, sticky='e')
        self.global_max_var= tk.StringVar(value="40")
        tk.Entry(fm_global, textvariable=self.global_max_var, width=5).grid(row=0, column=3, padx=5)

        fm_emp= tk.LabelFrame(root, text="Employees", padx=5, pady=5)
        fm_emp.pack(fill='x', padx=10, pady=5)

        tk.Label(fm_emp, text="# of employees:").grid(row=0, column=0, sticky='e')
        self.num_emp_var= tk.StringVar(value="3")
        tk.Entry(fm_emp, textvariable=self.num_emp_var, width=5).grid(row=0, column=1, padx=5)
        tk.Button(fm_emp, text="Generate Rows", command=self.gen_emp_rows).grid(row=0, column=2, padx=5)

        self.emp_rows_frame= tk.Frame(fm_emp)
        self.emp_rows_frame.grid(row=1, column=0, columnspan=3)
        self.emp_rows= []

        fm_forb= tk.LabelFrame(root, text="Forbidden Pairs", padx=5, pady=5)
        fm_forb.pack(fill='x', padx=10, pady=5)
        tk.Label(fm_forb, text="(emp1 emp2 per line)").pack(side='left')
        self.forbidden_text= tk.Text(fm_forb, width=40, height=5)
        self.forbidden_text.pack(side='left', padx=5)

        fm_bottom= tk.Frame(root)
        fm_bottom.pack(fill='both', expand=True, padx=10, pady=5)

        tk.Button(fm_bottom, text="Solve & Generate", command=self.solve_and_display).pack(side='left', padx=5)

        self.result_text= tk.Text(fm_bottom, width=60, height=10)
        self.result_text.pack(side='left', fill='both', expand=True, padx=5)

        self.export_btn= tk.Button(fm_bottom, text="Export CSV", command=self.export_csv, state='disabled')
        self.export_btn.pack(side='left', padx=5)

        self.schedule_data= None
        self.gen_emp_rows()

    def gen_emp_rows(self):
        for w in self.emp_rows_frame.winfo_children():
            w.destroy()
        self.emp_rows=[]
        try:
            n= int(self.num_emp_var.get())
        except:
            n=3
        for i in range(n):
            rowf= tk.Frame(self.emp_rows_frame)
            rowf.pack(fill='x', pady=3)
            name_var= tk.StringVar(value=f"Emp{i+1}")
            ds_var= tk.StringVar(value="2")

            tk.Label(rowf, text=f"Employee #{i+1}:").pack(side='left')
            tk.Label(rowf, text="Name").pack(side='left')
            e_name= tk.Entry(rowf, textvariable=name_var, width=10)
            e_name.pack(side='left', padx=5)

            tk.Label(rowf, text="Desired shift").pack(side='left')
            e_ds= tk.Entry(rowf, textvariable=ds_var, width=3)
            e_ds.pack(side='left', padx=5)

            tk.Label(rowf, text="Availability\n(day start end or 'all'; lines)").pack(side='left', padx=5)
            txt_avail= tk.Text(rowf, width=30, height=3)
            txt_avail.pack(side='left', padx=5)

            self.emp_rows.append((name_var, ds_var, txt_avail))

    def solve_and_display(self):
        self.result_text.delete("1.0","end")
        self.export_btn.config(state='disabled')
        try:
            gmin= int(self.global_min_var.get())
            gmax= int(self.global_max_var.get())
        except:
            messagebox.showerror("Error","Invalid global min/max hours.")
            return

        employees_availability={}
        desired_shifts_map={}

        for (name_var, ds_var, txt_avail) in self.emp_rows:
            name= name_var.get().strip()
            if not name:
                messagebox.showerror("Error", "Empty employee name.")
                return
            try:
                ds= int(ds_var.get())
            except:
                messagebox.showerror("Error", f"{name}: invalid desired shift.")
                return
            lines= txt_avail.get("1.0","end").strip().split("\n")
            avset= set()
            for ln in lines:
                ln= ln.strip().lower()
                if not ln or ln=="done":
                    continue
                if ln=="all":
                    avset= set(range(7*12))
                    break
                parts= ln.split()
                if len(parts)==3:
                    try:
                        day= int(parts[0])
                        sh= int(parts[1])
                        eh= int(parts[2])
                        if 1<=day<=7 and 11<=sh<eh<=23:
                            d_idx= day-1
                            for hour in range(sh,eh):
                                slot_id= d_idx*12 + (hour-11)
                                avset.add(slot_id)
                        else:
                            messagebox.showerror("Error", f"Invalid line '{ln}' for {name}. out of range.")
                            return
                    except:
                        messagebox.showerror("Error", f"Invalid line '{ln}' for {name}. expect 'day start end'")
                        return
                else:
                    messagebox.showerror("Error", f"Invalid line '{ln}' for {name}.")
                    return
            employees_availability[name]= avset
            desired_shifts_map[name]= ds

        lines_forb= self.forbidden_text.get("1.0","end").strip().split("\n")
        forbidden_pairs=[]
        for ln in lines_forb:
            ln= ln.strip()
            if not ln:
                continue
            parts= ln.split()
            if len(parts)==2:
                forbidden_pairs.append((parts[0], parts[1]))
            else:
                messagebox.showerror("Error", f"Invalid forbidden pair '{ln}' (emp1 emp2).")
                return

        # 1) Attempt main solve
        schedule_main, solver_main, status_main, tsv_main= solve_scheduling_main(
            employees_availability,
            global_min_hours=gmin,
            global_max_hours=gmax,
            forbidden_pairs= forbidden_pairs,
            coverage_needed=2,
            desired_shifts= desired_shifts_map
        )
        if schedule_main is not None:
            self.result_text.insert("end", f"Main solve found a perfect coverage=2 schedule! status={status_main}\n")
            for e in employees_availability.keys():
                hrs= len(schedule_main[e])
                got= solver_main.Value(tsv_main[e])
                self.result_text.insert("end", f"  {e}: {hrs} hours, desired={desired_shifts_map[e]}, got={got}\n")
            self.schedule_data= (schedule_main, solver_main, tsv_main, desired_shifts_map)
            self.export_btn.config(state='normal')
        else:
            self.result_text.insert("end", f"Main solve infeasible. status={status_main}\n")
            self.result_text.insert("end", "Attempting fallback partial coverage solve...\n")
            # 2) attempt fallback
            schedule_fb, solver_fb, status_fb, tsv_fb= solve_scheduling_fallback(
                employees_availability,
                global_min_hours=gmin,
                global_max_hours=gmax,
                forbidden_pairs= forbidden_pairs,
                desired_shifts= desired_shifts_map
            )
            if schedule_fb is None:
                self.result_text.insert("end", f"Fallback also infeasible. status={status_fb}\nNo schedule can be produced.\n")
            else:
                self.result_text.insert("end", f"Fallback partial coverage solve success! status={status_fb}\n")
                # print coverage
                # sum coverage
                coverage_sum=0
                for e in schedule_fb:
                    coverage_sum += len(schedule_fb[e])
                self.result_text.insert("end", f"Coverage used = {coverage_sum} (some hours may be partially or not covered)\n")
                for e in employees_availability.keys():
                    hrs= len(schedule_fb[e])
                    got= solver_fb.Value(tsv_fb[e])
                    self.result_text.insert("end", f"  {e}: {hrs} hours, desired={desired_shifts_map[e]}, got={got}\n")
                self.schedule_data= (schedule_fb, solver_fb, tsv_fb, desired_shifts_map)
                self.export_btn.config(state='normal')

    def export_csv(self):
        if not self.schedule_data:
            messagebox.showinfo("No data","No schedule to export.")
            return
        filename= filedialog.asksaveasfilename(title="Save CSV as...", defaultextension=".csv")
        if not filename:
            return
        schedule, solver, total_shifts_vars, desired_shifts_map= self.schedule_data
        export_schedule_to_csv(schedule, list(schedule.keys()), solver, total_shifts_vars, desired_shifts_map, filename)
        messagebox.showinfo("Export","CSV exported successfully!")


def main():
    root= tk.Tk()
    app= SinglePageSchedulerGUI(root)
    root.mainloop()

if __name__=="__main__":
    main()