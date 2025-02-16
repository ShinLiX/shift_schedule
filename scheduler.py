#!/usr/bin/env python3
import tkinter as tk
from tkinter import messagebox, filedialog
import csv
from ortools.sat.python import cp_model

# ========================== 1) Support Functions ==========================

def compute_total_weekly_availability(employees_availability):
    """Return {employee: total # of available weekly slots (0..83)}."""
    return { e: len(avail) for e, avail in employees_availability.items() }

def is_contiguous_ones(bits):
    """Check if bits is exactly one contiguous block of 1s, or all zero."""
    idxs = [i for i,b in enumerate(bits) if b==1]
    if not idxs:
        return True
    first = idxs[0]
    last  = idxs[-1]
    return (last - first + 1) == len(idxs)

def generate_day_patterns(availability_mask, min_block=6, max_block=8):
    """
    For a single day of 12 hours => sum=0 => off,
    if sum(availability_mask)<6 => fill smaller block if contiguous or off,
    else => one contiguous block in [6..8] or off.
    Return patterns of length=13 => 12 bits + 1 bit for 'workedBit' (did we work?).
    """
    n = sum(availability_mask)
    valid_pats = []
    for mask in range(1<<12):
        pat = [(mask >> i)&1 for i in range(12)]
        # Must be subset of availability
        if any(pat[i]==1 and availability_mask[i]==0 for i in range(12)):
            continue

        total_ones = sum(pat)
        worked_bit = 1 if total_ones>0 else 0

        if n < min_block:
            # if daily availability <6 => either sum=0 (off) or fill entire block= n if contiguous
            if total_ones == 0:
                valid_pats.append(pat+[0])
            elif total_ones == n:
                # must match availability + contiguous
                same = True
                for i in range(12):
                    if pat[i] != availability_mask[i]:
                        same=False
                        break
                if same and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
        else:
            if total_ones == 0:
                valid_pats.append(pat+[0])
            else:
                if (min_block<= total_ones <= max_block) and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
    return valid_pats


# ========================== 2) Main Solve: coverage=2 EXACT ==========================
def solve_scheduling_main(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    desired_shifts=None
):
    """
    EXACT coverage=2 approach. If infeasible => returns None.
    Single-shift logic, min/max hours, forbidden pairs, desired shifts etc.
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees_availability }

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    x={}
    for e in employees_availability:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage=2 => sum(x[e,s])=2 each slot
    for s in range(total_slots):
        model.Add( sum(x[(e,s)] for e in employees_availability) == 2 )

    # if not in availability => x[e,s]=0
    for e,availset in employees_availability.items():
        for s in range(total_slots):
            if s not in availset:
                model.Add(x[(e,s)]==0)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees_availability:
        ds= desired_shifts[e]
        eff_min= 0 if (tot_avail[e]<global_min_hours or ds*8<global_min_hours) else global_min_hours
        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add(weekly_hours[e] >= eff_min)
        model.Add(weekly_hours[e] <= global_max_hours)

    # forbidden pairs => x[e1,s]+ x[e2,s]<=1
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees_availability and emp2 in employees_availability:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)] <= 1)

    # day-based single shift pattern
    day_pattern_index={}
    day_worked={}
    all_patterns={}
    for e in employees_availability:
        all_patterns[e]={}
        for d in range(days):
            avmask=[]
            startSlot= d*slots_per_day
            for i in range(slots_per_day):
                slot_id= startSlot + i
                avmask.append(1 if slot_id in employees_availability[e] else 0)

            valid_pats= generate_day_patterns(avmask,6,8)
            patIdx= model.NewIntVar(0,len(valid_pats)-1, f"patIndex_{e}_{d}")
            day_pattern_index[(e,d)] = patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    for e in employees_availability:
        for d in range(days):
            day_vars= [ x[(e,d*12 + i)] for i in range(12)]
            plus= day_vars + [ day_worked[(e,d)], day_pattern_index[(e,d)] ]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12] + [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # total_shifts => sum(day_worked[e,d])
    total_shifts={}
    diffs=[]
    for e in employees_availability:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"total_shifts_{e}")
        model.Add( sum(sumDays) == total_shifts[e])

    # handle desired shifts
    for e in employees_availability:
        ds= desired_shifts[e]
        if ds in [1,2]:
            model.Add( total_shifts[e]== ds )
        else:
            model.Add( total_shifts[e]<= ds )
            diffUp= model.NewIntVar(0,7,f"diffUp_{e}")
            diffDown= model.NewIntVar(0,7,f"diffDown_{e}")
            model.Add( total_shifts[e] - ds <= diffUp )
            model.Add( ds - total_shifts[e] <= diffDown )
            diffs.append(diffUp)
            diffs.append(diffDown)

    model.Minimize( sum(diffs) )

    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    # build schedule
    schedule= { e:[] for e in employees_availability }
    for e in employees_availability:
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
    If main solve fails, do partial coverage => sum(x[e,s])<=2, maximize total coverage. 
    Keep single-shift constraints, min/max hours, forbidden pairs, desired shifts, etc.
    Some hours might remain uncovered => occupant=NA in CSV
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees_availability }

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    x={}
    for e in employees_availability:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage <=2
    for s in range(total_slots):
        model.Add( sum(x[(e,s)] for e in employees_availability ) <= 2 )

    # availability
    for e,avset in employees_availability.items():
        for s in range(total_slots):
            if s not in avset:
                model.Add(x[(e,s)]==0)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees_availability:
        ds= desired_shifts[e]
        eff_min= 0 if (tot_avail[e]<global_min_hours or ds*8<global_min_hours) else global_min_hours
        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add(weekly_hours[e]>= eff_min)
        model.Add(weekly_hours[e]<= global_max_hours)

    # forbidden
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees_availability and emp2 in employees_availability:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)]<=1)

    # day-based single shift
    day_pattern_index={}
    day_worked={}
    all_patterns={}
    for e in employees_availability:
        all_patterns[e]={}
        for d in range(days):
            avmask=[]
            startSlot= d*12
            for i in range(12):
                slot_id= startSlot + i
                avmask.append(1 if slot_id in employees_availability[e] else 0)
            valid_pats= generate_day_patterns(avmask,6,8)
            patIdx= model.NewIntVar(0,len(valid_pats)-1,f"fb_patIndex_{e}_{d}")
            day_pattern_index[(e,d)]= patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"fb_dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    for e in employees_availability:
        for d in range(days):
            day_vars= [x[(e,d*12 + i)] for i in range(12)]
            plus= day_vars + [ day_worked[(e,d)], day_pattern_index[(e,d)] ]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12]+ [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # desired shifts => min difference
    total_shifts={}
    diffs=[]
    for e in employees_availability:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"fb_total_shifts_{e}")
        model.Add( sum(sumDays)== total_shifts[e])
        ds= desired_shifts[e]
        if ds in [1,2]:
            model.Add( total_shifts[e]== ds )
        else:
            model.Add( total_shifts[e]<= ds )
            diffUp= model.NewIntVar(0,7,f"fb_diffUp_{e}")
            diffDown= model.NewIntVar(0,7,f"fb_diffDown_{e}")
            model.Add( total_shifts[e]- ds <= diffUp )
            model.Add( ds - total_shifts[e] <= diffDown )
            diffs.append(diffUp)
            diffs.append(diffDown)

    # define coverage_sum => sum of all x[e,s]
    coverage_sum= model.NewIntVar(0,7*12*2,"coverage_sum")
    model.Add( coverage_sum == sum(x[(e,s)] for e in employees_availability for s in range(total_slots)) )

    # Weighted approach => maximize coverage_sum*1000 - sum(diffs) 
    # => prioritize coverage first, then shift difference
    model.Maximize( coverage_sum*1000 - sum(diffs) )

    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    schedule= { e:[] for e in employees_availability}
    for e in employees_availability:
        for s in range(total_slots):
            if solver.Value(x[(e,s)])==1:
                schedule[e].append(s)
    return schedule, solver, status, total_shifts


# ========================== 4) CSV Export w/ NA for partial coverage ==========================
def export_schedule_to_csv(schedule, employees, solver, total_shifts_vars, desired_shifts, filename="schedule.csv"):
    days=7
    slots_per_day=12
    day_names= ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]

    # build occupant info
    day_hour_emps= { d: {h: [] for h in range(slots_per_day)} for d in range(days)}
    coverage_count= { (d,h):0 for d in range(days) for h in range(slots_per_day)}

    for e in employees:
        for slot_id in schedule[e]:
            d= slot_id// slots_per_day
            h= slot_id% slots_per_day
            day_hour_emps[d][h].append(e)

    # occupant continuity
    day_col_occ= { d:{1:["" for _ in range(slots_per_day)],
                      2:["" for _ in range(slots_per_day)]}
                   for d in range(days)}

    for d in range(days):
        occupant1=None
        occupant2=None
        for h in range(slots_per_day):
            assigned= sorted(day_hour_emps[d][h])
            coverage_count[(d,h)]= len(assigned)  # could be 0..2
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

    # build occupant blocks so we only print occupant once
    day_col_blocks= { d:{1:[], 2:[]} for d in range(days)}
    for d in range(days):
        for col in [1,2]:
            occupant_list= day_col_occ[d][col]
            h=0
            while h< slots_per_day:
                name= occupant_list[h]
                if name=="":
                    day_col_blocks[d][col].append((h,h,""))
                    h+=1
                else:
                    startH=h
                    while (h+1<slots_per_day) and occupant_list[h+1]== name:
                        h+=1
                    endH=h
                    day_col_blocks[d][col].append((startH,endH,name))
                    h+=1

    occupantPrint= { d:{col:[""]*slots_per_day for col in [1,2]} for d in range(days)}
    for d in range(days):
        for col in [1,2]:
            for (startH,endH,name) in day_col_blocks[d][col]:
                occupantPrint[d][col][startH]= name

    # build table
    headers=["Time"]
    for d in range(days):
        headers.append(f"{day_names[d]}-1")
        headers.append(f"{day_names[d]}-2")

    rows=[]
    for hour in range(slots_per_day):
        timeLabel= f"{11+hour}:00 - {12+hour}:00"
        row_cells= [timeLabel]
        for d in range(days):
            cov= coverage_count[(d,hour)]
            c1= occupantPrint[d][1][hour]
            c2= occupantPrint[d][2][hour]
            if cov==0:
                c1="NA"
                c2="NA"
            elif cov==1:
                # occupant2 => NA
                if not c1 and c2:
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

# ========================== 5) Single-Page Tkinter GUI with fallback & partial coverage + "NA" + coverage=2 ==========================
class SinglePageSchedulerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Employee Scheduler - Single Page with fallback coverage")

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

            tk.Label(rowf, text="Availability\n(day start end or 'all'; lines)\n(We do half-open => end - 1 if needed)").pack(side='left', padx=5)
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

                        # subtract 1 from end => half-open range [sh..eh-1)
                        eh_adjusted= eh - 1
                        if eh_adjusted < sh:
                            # no interval
                            continue

                        if 1<= day<=7:
                            d_idx= day-1
                            # clamp to [11..22] if needed
                            if sh<11: 
                                sh= 11
                            if eh_adjusted> 22: 
                                eh_adjusted=22
                            for hour in range(sh, eh_adjusted+1):
                                slot_index= hour - 11
                                if 0 <= slot_index < 12:
                                    slot_id= d_idx*12 + slot_index
                                    avset.add(slot_id)
                        else:
                            messagebox.showerror("Error", f"Invalid line '{ln}' => day out of range 1..7.")
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

        # 1) Attempt main coverage=2 solve
        schedule_main, solver_main, status_main, tsv_main= solve_scheduling_main(
            employees_availability,
            global_min_hours=gmin,
            global_max_hours=gmax,
            forbidden_pairs= forbidden_pairs,
            desired_shifts= desired_shifts_map
        )
        if schedule_main is not None:
            self.result_text.insert("end", f"Main solve found coverage=2 schedule! status={status_main}\n")
            for e in employees_availability:
                hrs= len(schedule_main[e])
                got= solver_main.Value(tsv_main[e])
                self.result_text.insert("end", f"  {e}: {hrs} hours, desired={desired_shifts_map[e]}, got={got}\n")
            self.schedule_data= (schedule_main, solver_main, tsv_main, desired_shifts_map)
            self.export_btn.config(state='normal')
        else:
            self.result_text.insert("end", f"Main solve infeasible. status={status_main}\n")
            self.result_text.insert("end", "Attempting fallback partial coverage solve...\n")
            # 2) fallback coverage <=2, maximize coverage
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
                self.result_text.insert("end", f"Fallback partial coverage success! status={status_fb}\n")
                coverage_sum= sum(len(schedule_fb[e]) for e in employees_availability)
                self.result_text.insert("end", f"Coverage used = {coverage_sum}\n(some hours partially or not covered => NA)\n")
                for e in employees_availability:
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