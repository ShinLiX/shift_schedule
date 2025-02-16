#!/usr/bin/env python3
import tkinter as tk
from tkinter import messagebox, filedialog
import csv
from ortools.sat.python import cp_model

# ========================== 1) Support / Common ==========================

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
    For a single day (12 hours => 0..11):
      - if sum(availability_mask)<6 => fill entire smaller block if contiguous or off
      - else => one contiguous block in [6..8] or off
      - return patterns of length=13 => 12 bits + 1 workedBit
    """
    n = sum(availability_mask)
    valid_pats = []
    for mask in range(1<<12):
        pat = [(mask >> i)&1 for i in range(12)]
        # Must be subset
        if any(pat[i]==1 and availability_mask[i]==0 for i in range(12)):
            continue

        total_ones = sum(pat)
        worked_bit = 1 if total_ones>0 else 0

        if n< min_block:
            if total_ones==0:
                valid_pats.append(pat+[0])
            elif total_ones==n:
                # must match availability exactly & be contiguous
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
                if (min_block<= total_ones <= max_block) and is_contiguous_ones(pat):
                    valid_pats.append(pat+[1])
    return valid_pats

# ========================== 2) Main Solve => coverage=2 EXACT ==========================
def solve_scheduling_main(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    desired_shifts=None
):
    """
    EXACT coverage=2 each hour. If infeasible => returns None.
    Single-shift logic, min/max hours, forbidden pairs, desired_shifts in [1,2 => forced, else <=].
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees_availability }

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    # x[e,s]=1 => e works slot s
    x={}
    for e in employees_availability:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # coverage=2 => sum(x[e,s])=2 each slot
    for s in range(total_slots):
        model.Add( sum(x[(e,s)] for e in employees_availability) == 2 )

    # if not in availability => x[e,s]=0
    for e,avset in employees_availability.items():
        for s in range(total_slots):
            if s not in avset:
                model.Add(x[(e,s)]==0)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees_availability:
        ds= desired_shifts[e]
        eff_min= 0
        if not (tot_avail[e]< global_min_hours or ds*8< global_min_hours):
            eff_min= global_min_hours

        weekly_hours[e]= sum(x[(e,s)] for s in range(total_slots))
        model.Add( weekly_hours[e]>= eff_min )
        model.Add( weekly_hours[e]<= global_max_hours)

    # forbidden pairs => x[e1,s]+ x[e2,s] <=1
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees_availability and emp2 in employees_availability:
            for s in range(total_slots):
                model.Add(x[(emp1,s)] + x[(emp2,s)]<=1)

    # day-based single shift => generate patterns
    day_pattern_index={}
    day_worked={}
    all_patterns={}
    for e in employees_availability:
        all_patterns[e]={}
        for d in range(days):
            base= d*12
            avmask=[]
            for i in range(12):
                slot_id= base + i
                avmask.append(1 if slot_id in employees_availability[e] else 0)

            valid_pats= generate_day_patterns(avmask,6,8)
            patIdx= model.NewIntVar(0, len(valid_pats)-1, f"patIndex_{e}_{d}")
            day_pattern_index[(e,d)] = patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    for e in employees_availability:
        for d in range(days):
            day_vars= [ x[(e,d*12 + i)] for i in range(12) ]
            plus= day_vars + [day_worked[(e,d)], day_pattern_index[(e,d)]]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12] + [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # total_shifts => sum day_worked
    total_shifts={}
    diffs=[]
    for e in employees_availability:
        sumDays= [ day_worked[(e,d)] for d in range(days)]
        total_shifts[e]= model.NewIntVar(0,7,f"total_shifts_{e}")
        model.Add( sum(sumDays)== total_shifts[e])

    # desired_shifts => forced if in [1,2], else <= ds => min difference
    for e in employees_availability:
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
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    schedule= { e:[] for e in employees_availability}
    for e in employees_availability:
        for s in range(total_slots):
            if solver.Value(x[(e,s)])==1:
                schedule[e].append(s)
    return schedule, solver, status, total_shifts


# ========================== 3) Fallback: coverage ≤2, single contiguous block or full day (12) per day ==========================
def solve_scheduling_fallback(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    desired_shifts=None
):
    """
    If main solve fails, fallback coverage <=2. 
    Each day => coverage_mask in [0..1]^12 with sum=0, sum=12, or single contiguous block => coverage=2 in that block, else 0.
    Also keep single-shift logic, min/max hours, forbidden pairs, desired shifts in [1,2 => forced, else <=].
    Maximize coverage minus shift diffs.
    """
    if forbidden_pairs is None:
        forbidden_pairs=[]
    if desired_shifts is None:
        desired_shifts={ e:0 for e in employees_availability }

    model= cp_model.CpModel()
    days=7
    slots_per_day=12
    total_slots= days*slots_per_day

    # coverage_mask[(d,h)] in {0,1} => if 1 => coverage=2 that hour, else coverage=0
    coverage_mask={}
    for d in range(days):
        for h in range(slots_per_day):
            coverage_mask[(d,h)] = model.NewBoolVar(f"covMask_{d}_{h}")

    # day-based => sum(coverage_mask[d,h]) in {0,12 or contiguous block}
    for d in range(days):
        day_mask_vars= [ coverage_mask[(d,h)] for h in range(slots_per_day)]
        valid_daycov=[]
        # enumerate all 2^12 => if sum=0 => all off, sum=12 => all on, or sum>0<12 => must be contiguous
        for mask in range(1<<12):
            pat= [(mask >> i)&1 for i in range(slots_per_day)]
            s_cov= sum(pat)
            if s_cov==0 or s_cov==12:
                valid_daycov.append(pat)
            else:
                # must be contiguous if sum>0 <12
                if is_contiguous_ones(pat):
                    valid_daycov.append(pat)
        model.AddAllowedAssignments(day_mask_vars, valid_daycov)

    # define x[e,s] => employee e works slot s => if coverage_mask[d,h]=1 => sum(x[e,d*12+h])=2, else=0
    x={}
    for e in employees_availability:
        for s in range(total_slots):
            x[(e,s)] = model.NewBoolVar(f'x_{e}_{s}')

    # link coverage_mask => sum(x[e,s])= coverage_mask*(2)
    for d in range(days):
        for h in range(slots_per_day):
            s_id= d*slots_per_day + h
            # sum(x[e,s_id])= coverage_mask[d,h]*2
            model.Add( sum(x[(e,s_id)] for e in employees_availability) == 2* coverage_mask[(d,h)] )

    # if not in availability => x[e,s]=0
    for e,avset in employees_availability.items():
        for s in range(total_slots):
            if s not in avset:
                model.Add(x[(e,s)]==0)

    # single-shift logic
    day_pattern_index={}
    day_worked={}
    all_patterns={}
    for e in employees_availability:
        all_patterns[e]={}
        for d in range(days):
            base= d*slots_per_day
            avmask=[]
            for i in range(slots_per_day):
                slot_id= base + i
                avmask.append(1 if slot_id in employees_availability[e] else 0)
            valid_pats= generate_day_patterns(avmask,6,8)
            patIdx= model.NewIntVar(0,len(valid_pats)-1, f"fbEmp_patIndex_{e}_{d}")
            day_pattern_index[(e,d)] = patIdx
            day_worked[(e,d)] = model.NewBoolVar(f"fbEmp_dayWorked_{e}_{d}")
            all_patterns[e][d]= valid_pats

    for e in employees_availability:
        for d in range(days):
            day_vars_emp= [ x[(e,d*12 + i)] for i in range(12)]
            plus= day_vars_emp + [ day_worked[(e,d)], day_pattern_index[(e,d)] ]
            table=[]
            for idx, pat in enumerate(all_patterns[e][d]):
                row= pat[:12] + [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # min/max hours
    tot_avail= compute_total_weekly_availability(employees_availability)
    weekly_hours={}
    for e in employees_availability:
        ds= desired_shifts[e]
        eff_min=0
        if not (tot_avail[e]<global_min_hours or ds*8<global_min_hours):
            eff_min= global_min_hours
        wh= sum( x[(e,s)] for s in range(total_slots))
        weekly_hours[e]= wh
        model.Add( wh>= eff_min )
        model.Add( wh<= global_max_hours)

    # forbidden pairs => x[e1,s]+ x[e2,s]<=1
    for (emp1,emp2) in forbidden_pairs:
        if emp1 in employees_availability and emp2 in employees_availability:
            for s in range(total_slots):
                model.Add( x[(emp1,s)] + x[(emp2,s)]<=1 )

    # desired shifts => sum day_worked => min difference
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
            model.Add( ds- total_shifts[e] <= diffDown )
            diffs.append(diffUp)
            diffs.append(diffDown)

    # coverage_sum => sum( coverage_mask*(2) ) => we want to maximize coverage
    coverage_sum= model.NewIntVar(0,7*12*2,"coverage_sum")
    model.Add( coverage_sum == sum( coverage_mask[(d,h)]*2 for d in range(days) for h in range(slots_per_day)) )

    # Weighted approach => coverage_sum*1000 - sum(diffs)
    model.Maximize( coverage_sum*1000 - sum(diffs))

    solver= cp_model.CpSolver()
    status= solver.Solve(model)
    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        return None, None, status, None

    # build schedule => x[e,s]==1 => e works slot s
    schedule= { e:[] for e in employees_availability}
    for e in employees_availability:
        for s in range(total_slots):
            if solver.Value( x[(e,s)])==1:
                schedule[e].append(s)
    return schedule, solver, status, total_shifts


# ========================== 4) CSV export => occupant repeated each hour, "NA" for partial coverage ==========================
def export_schedule_to_csv(
    schedule, 
    employees, 
    solver, 
    total_shifts_vars, 
    desired_shifts, 
    filename="schedule.csv"
):
    """
    We skip occupant continuity merges. For each day/hour => occupant1, occupant2 => if coverage=0 => both=NA, if coverage=1 => occupant2=NA
    leftover NA hours are contiguous in fallback, occupant repeated each hour
    """
    days=7
    slots_per_day=12
    day_names= ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]

    day_hour_emps= { d:{h: [] for h in range(slots_per_day)} for d in range(days)}
    for e in employees:
        for s in schedule[e]:
            d= s// slots_per_day
            h= s%  slots_per_day
            day_hour_emps[d][h].append(e)

    occupant1={}
    occupant2={}
    for d in range(days):
        for h in range(slots_per_day):
            assigned= sorted(day_hour_emps[d][h])
            if len(assigned)==0:
                occupant1[(d,h)]="NA"
                occupant2[(d,h)]="NA"
            elif len(assigned)==1:
                occupant1[(d,h)]= assigned[0]
                occupant2[(d,h)]="NA"
            else:
                occupant1[(d,h)]= assigned[0]
                occupant2[(d,h)]= assigned[1]

    headers=["Time"]
    for d in range(days):
        headers.append(f"{day_names[d]}-1")
        headers.append(f"{day_names[d]}-2")

    rows=[]
    for hour in range(slots_per_day):
        time_str= f"{11+hour}:00 - {12+hour}:00"
        row_cells=[time_str]
        for d in range(days):
            row_cells.append( occupant1[(d,hour)] )
            row_cells.append( occupant2[(d,hour)] )
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
        for row in rows:
            writer.writerow(row)

    print(f"Exported schedule to '{filename}'.")


# ========================== 5) Single-Page Tkinter GUI with fallback coverage ==========================
class SinglePageSchedulerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Scheduler with leftover NA grouped, occupant repeated each hour")

        # Global constraints
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
        self.emp_rows=[]

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

            tk.Label(rowf, text="Availability\n(day start end)\nWe do half-open => end-1").pack(side='left', padx=5)
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

        # parse employees
        for (name_var, ds_var, txt_avail) in self.emp_rows:
            name= name_var.get().strip()
            if not name:
                messagebox.showerror("Error","Empty employee name.")
                return
            try:
                ds= int(ds_var.get())
            except:
                messagebox.showerror("Error",f"{name}: invalid desired shift.")
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

                        # half-open => subtract 1 from end
                        eh_adj= eh -1
                        if eh_adj< sh:
                            continue
                        if day<1 or day>7:
                            messagebox.showerror("Error", f"Invalid day in line '{ln}' => out of range [1..7]")
                            return
                        d_idx= day-1
                        # clamp to [11..22]
                        if sh<11: sh=11
                        if eh_adj>22: eh_adj=22
                        for hour in range(sh, eh_adj+1):
                            slot_index= hour-11
                            if 0<= slot_index<12:
                                slot_id= d_idx*12 + slot_index
                                avset.add(slot_id)
                    except:
                        messagebox.showerror("Error", f"Invalid line '{ln}' for {name}. Expect 'day start end'.")
                        return
                else:
                    messagebox.showerror("Error", f"Invalid line '{ln}' for {name}.")
                    return
            employees_availability[name]= avset
            desired_shifts_map[name]= ds

        # parse forbidden
        forb_lines= self.forbidden_text.get("1.0","end").strip().split("\n")
        forbidden_pairs=[]
        for ln in forb_lines:
            ln= ln.strip()
            if not ln:
                continue
            pair= ln.split()
            if len(pair)==2:
                forbidden_pairs.append((pair[0], pair[1]))
            else:
                messagebox.showerror("Error", f"Invalid forbidden pair '{ln}'.")
                return

        # 1) main solve coverage=2
        schedule_main, solver_main, status_main, tsv_main= solve_scheduling_main(
            employees_availability,
            global_min_hours=gmin,
            global_max_hours=gmax,
            forbidden_pairs= forbidden_pairs,
            desired_shifts= desired_shifts_map
        )
        if schedule_main is not None:
            self.result_text.insert("end", f"Main solve => coverage=2 success! status={status_main}\n")
            for e in employees_availability:
                hrs= len(schedule_main[e])
                got= solver_main.Value(tsv_main[e])
                self.result_text.insert("end", f"  {e}: {hrs} hrs, desired={desired_shifts_map[e]}, got={got}\n")
            self.schedule_data= (schedule_main, solver_main, tsv_main, desired_shifts_map)
            self.export_btn.config(state='normal')
        else:
            # fallback coverage => single contiguous block or full day => leftover NA grouped
            self.result_text.insert("end", f"Main solve infeasible => status={status_main}\n")
            self.result_text.insert("end", "Attempt fallback coverage...\n")
            schedule_fb, solver_fb, status_fb, tsv_fb= solve_scheduling_fallback(
                employees_availability,
                global_min_hours=gmin,
                global_max_hours=gmax,
                forbidden_pairs= forbidden_pairs,
                desired_shifts= desired_shifts_map
            )
            if schedule_fb is None:
                self.result_text.insert("end", f"Fallback also infeasible => status={status_fb}\nNo schedule.\n")
            else:
                self.result_text.insert("end", f"Fallback partial coverage success => status={status_fb}\n(Leftover NA hours grouped morning or night)\n")
                coverage_sum= sum(len(schedule_fb[e]) for e in employees_availability)
                self.result_text.insert("end", f"Coverage used= {coverage_sum}\n")
                for e in employees_availability:
                    hrs= len(schedule_fb[e])
                    got= solver_fb.Value(tsv_fb[e])
                    self.result_text.insert("end", f"  {e}: {hrs} hrs, desired={desired_shifts_map[e]}, got={got}\n")
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
        export_schedule_to_csv(
            schedule, 
            list(schedule.keys()), 
            solver, 
            total_shifts_vars, 
            desired_shifts_map, 
            filename
        )
        messagebox.showinfo("Export","CSV exported successfully!")


def main():
    root= tk.Tk()
    app= SinglePageSchedulerGUI(root)
    root.mainloop()

if __name__=="__main__":
    main()