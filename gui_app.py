# gui_app.py
"""
Single-Page Tkinter GUI with a 2-phase solve and fallback partial coverage.
"""

import tkinter as tk
from tkinter import messagebox, filedialog

from main_solver import solve_scheduling_main
from fallback_solver import solve_scheduling_fallback
from exporter import export_schedule_to_csv

class SinglePageSchedulerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Employee Scheduler - Single Page (Fallback Partial Coverage)")

        # 1) Global constraints frame
        fm_global = tk.LabelFrame(root, text="Global Constraints", padx=5, pady=5)
        fm_global.pack(fill='x', padx=10, pady=5)

        tk.Label(fm_global, text="Global min hours:").grid(row=0, column=0, sticky='e')
        self.global_min_var = tk.StringVar(value="16")
        tk.Entry(fm_global, textvariable=self.global_min_var, width=5).grid(row=0, column=1, padx=5)

        tk.Label(fm_global, text="Global max hours:").grid(row=0, column=2, sticky='e')
        self.global_max_var = tk.StringVar(value="40")
        tk.Entry(fm_global, textvariable=self.global_max_var, width=5).grid(row=0, column=3, padx=5)

        # 2) Employee frames
        fm_emp = tk.LabelFrame(root, text="Employees", padx=5, pady=5)
        fm_emp.pack(fill='x', padx=10, pady=5)

        tk.Label(fm_emp, text="# of employees:").grid(row=0, column=0, sticky='e')
        self.num_emp_var = tk.StringVar(value="3")
        tk.Entry(fm_emp, textvariable=self.num_emp_var, width=5).grid(row=0, column=1, padx=5)
        tk.Button(fm_emp, text="Generate Rows", command=self.gen_emp_rows).grid(row=0, column=2, padx=5)

        self.emp_rows_frame = tk.Frame(fm_emp)
        self.emp_rows_frame.grid(row=1, column=0, columnspan=3)
        self.emp_rows = []

        # 3) Forbidden pairs
        fm_forb = tk.LabelFrame(root, text="Forbidden Pairs", padx=5, pady=5)
        fm_forb.pack(fill='x', padx=10, pady=5)
        tk.Label(fm_forb, text="(emp1 emp2 per line)").pack(side='left')
        self.forbidden_text = tk.Text(fm_forb, width=40, height=5)
        self.forbidden_text.pack(side='left', padx=5)

        # 4) Solve & result area
        fm_bottom = tk.Frame(root)
        fm_bottom.pack(fill='both', expand=True, padx=10, pady=5)

        tk.Button(fm_bottom, text="Solve & Generate", command=self.solve_and_display).pack(side='left', padx=5)

        self.result_text = tk.Text(fm_bottom, width=60, height=10)
        self.result_text.pack(side='left', fill='both', expand=True, padx=5)

        self.export_btn = tk.Button(fm_bottom, text="Export CSV", command=self.export_csv, state='disabled')
        self.export_btn.pack(side='left', padx=5)

        self.schedule_data = None  # store final solve result
        self.gen_emp_rows()

    def gen_emp_rows(self):
        """Generate n employee lines based on self.num_emp_var."""
        for w in self.emp_rows_frame.winfo_children():
            w.destroy()
        self.emp_rows = []
        try:
            n = int(self.num_emp_var.get())
        except ValueError:
            n = 3

        for i in range(n):
            rowf = tk.Frame(self.emp_rows_frame)
            rowf.pack(fill='x', pady=3)

            name_var = tk.StringVar(value=f"Emp{i+1}")
            ds_var = tk.StringVar(value="2")

            tk.Label(rowf, text=f"Employee #{i+1}:").pack(side='left')
            tk.Label(rowf, text="Name").pack(side='left')
            e_name = tk.Entry(rowf, textvariable=name_var, width=10)
            e_name.pack(side='left', padx=5)

            tk.Label(rowf, text="Desired shift").pack(side='left')
            e_ds = tk.Entry(rowf, textvariable=ds_var, width=3)
            e_ds.pack(side='left', padx=5)

            tk.Label(rowf, text="Availability\n(day start end)\nUse 'all' or 'done' if needed").pack(side='left', padx=5)
            txt_avail = tk.Text(rowf, width=30, height=3)
            txt_avail.pack(side='left', padx=5)

            self.emp_rows.append((name_var, ds_var, txt_avail))

    def solve_and_display(self):
        """Handle the solve button click, run main solver then fallback if needed."""
        self.result_text.delete("1.0", "end")
        self.export_btn.config(state='disabled')

        # Parse global min/max
        try:
            gmin = int(self.global_min_var.get())
            gmax = int(self.global_max_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid global min/max hours.")
            return

        employees_availability = {}
        desired_shifts_map = {}

        # Parse employees
        for (name_var, ds_var, txt_avail) in self.emp_rows:
            name = name_var.get().strip()
            if not name:
                messagebox.showerror("Error", "Empty employee name.")
                return
            try:
                ds = int(ds_var.get())
            except ValueError:
                messagebox.showerror("Error", f"{name}: invalid desired shift.")
                return

            lines = txt_avail.get("1.0", "end").strip().split("\n")
            avset = set()
            for ln in lines:
                ln = ln.strip().lower()
                if not ln or ln == "done":
                    continue
                if ln == "all":
                    avset = set(range(7 * 12))
                    break
                parts = ln.split()
                if len(parts) == 3:
                    try:
                        day = int(parts[0])
                        sh = int(parts[1])
                        eh = int(parts[2])
                        # half-open => end = eh - 1
                        eh_adjusted = eh - 1
                        if eh_adjusted < sh:
                            continue
                        if day < 1 or day > 7:
                            messagebox.showerror("Error", f"Invalid day in line '{ln}' => day out of [1..7]")
                            return
                        d_idx = day - 1

                        # clamp to 11..22 if needed
                        if sh < 11:
                            sh = 11
                        if eh_adjusted > 22:
                            eh_adjusted = 22

                        for hour in range(sh, eh_adjusted + 1):
                            slot_index = hour - 11
                            if 0 <= slot_index < 12:
                                slot_id = d_idx * 12 + slot_index
                                avset.add(slot_id)
                    except ValueError:
                        messagebox.showerror("Error", f"Invalid line '{ln}' for {name}. Expect 'day start end'.")
                        return
                else:
                    messagebox.showerror("Error", f"Invalid line '{ln}' for {name}.")
                    return

            employees_availability[name] = avset
            desired_shifts_map[name] = ds

        # Parse forbidden pairs
        forb_lines = self.forbidden_text.get("1.0", "end").strip().split("\n")
        forbidden_pairs = []
        for ln in forb_lines:
            ln = ln.strip()
            if not ln:
                continue
            pair = ln.split()
            if len(pair) == 2:
                forbidden_pairs.append((pair[0], pair[1]))
            else:
                messagebox.showerror("Error", f"Invalid forbidden pair '{ln}'. Expect 'emp1 emp2'.")
                return

        # 1) Main solve => coverage=2
        schedule_main, solver_main, status_main, tsv_main = solve_scheduling_main(
            employees_availability,
            global_min_hours=gmin,
            global_max_hours=gmax,
            forbidden_pairs=forbidden_pairs,
            desired_shifts=desired_shifts_map
        )

        if schedule_main is not None:
            # success coverage=2
            self.result_text.insert("end", f"Main solve => coverage=2 success! status={status_main}\n")
            for e in employees_availability:
                hrs = len(schedule_main[e])
                got = solver_main.Value(tsv_main[e])
                self.result_text.insert("end", f"  {e}: {hrs} hours, desired={desired_shifts_map[e]}, got={got}\n")
            self.schedule_data = (schedule_main, solver_main, tsv_main, desired_shifts_map)
            self.export_btn.config(state='normal')
        else:
            # fallback partial coverage
            self.result_text.insert("end", f"Main solve infeasible. status={status_main}\n")
            self.result_text.insert("end", "Attempting fallback partial coverage solve...\n")
            schedule_fb, solver_fb, status_fb, tsv_fb = solve_scheduling_fallback(
                employees_availability,
                global_min_hours=gmin,
                global_max_hours=gmax,
                forbidden_pairs=forbidden_pairs,
                desired_shifts=desired_shifts_map
            )
            if schedule_fb is None:
                self.result_text.insert("end", f"Fallback also infeasible => status={status_fb}\nNo schedule produced.\n")
            else:
                self.result_text.insert("end", f"Fallback partial coverage success => status={status_fb}\n")
                coverage_sum = sum(len(schedule_fb[e]) for e in employees_availability)
                self.result_text.insert("end", f"Coverage used={coverage_sum}. Some hours partial => 'NA' in CSV.\n")
                for e in employees_availability:
                    hrs = len(schedule_fb[e])
                    got = solver_fb.Value(tsv_fb[e])
                    self.result_text.insert("end", f"  {e}: {hrs} hours, desired={desired_shifts_map[e]}, got={got}\n")
                self.schedule_data = (schedule_fb, solver_fb, tsv_fb, desired_shifts_map)
                self.export_btn.config(state='normal')

    def export_csv(self):
        """Export the schedule data to a CSV file."""
        if not self.schedule_data:
            messagebox.showinfo("No data", "No schedule to export.")
            return

        filename = filedialog.asksaveasfilename(
            title="Save CSV as...",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not filename:
            return

        schedule, solver, total_shifts_vars, desired_shifts_map = self.schedule_data
        export_schedule_to_csv(
            schedule,
            list(schedule.keys()),
            solver,
            total_shifts_vars,
            desired_shifts_map,
            filename
        )
        messagebox.showinfo("Export", "CSV exported successfully!")