# main_solver.py
"""
Main solver for the scheduling problem with EXACT coverage=2.
"""

from ortools.sat.python import cp_model
from support import (
    compute_total_weekly_availability,
    generate_day_patterns
)

def solve_scheduling_main(
    employees_availability,
    global_min_hours=16,
    global_max_hours=40,
    forbidden_pairs=None,
    desired_shifts=None
):
    """
    Solve the scheduling problem requiring EXACT coverage=2 for each hour.
    
    If infeasible => returns (None, None, status, None).
    Otherwise => (schedule, solver, status, total_shifts_var_map).
    """
    if forbidden_pairs is None:
        forbidden_pairs = []
    if desired_shifts is None:
        desired_shifts = {e: 0 for e in employees_availability}

    model = cp_model.CpModel()

    days = 7
    slots_per_day = 12
    total_slots = days * slots_per_day

    # Decision variables: x[e,s] = 1 => employee e works slot s
    x = {}
    for e in employees_availability:
        for s in range(total_slots):
            x[(e, s)] = model.NewBoolVar(f'x_{e}_{s}')

    # Coverage=2 => sum(x[e,s]) == 2 for each slot
    for s in range(total_slots):
        model.Add(sum(x[(e, s)] for e in employees_availability) == 2)

    # If not available => x[e,s] = 0
    for e, avset in employees_availability.items():
        for s in range(total_slots):
            if s not in avset:
                model.Add(x[(e, s)] == 0)

    # Min/Max hours
    tot_avail = compute_total_weekly_availability(employees_availability)
    weekly_hours = {}
    for e in employees_availability:
        ds = desired_shifts[e]
        eff_min = 0
        # Check if we should enforce global_min_hours
        if not (tot_avail[e] < global_min_hours or ds * 8 < global_min_hours):
            eff_min = global_min_hours

        weekly_hours[e] = sum(x[(e, s)] for s in range(total_slots))
        model.Add(weekly_hours[e] >= eff_min)
        model.Add(weekly_hours[e] <= global_max_hours)

    # Forbidden pairs => can't share the same slot
    for (emp1, emp2) in forbidden_pairs:
        if emp1 in employees_availability and emp2 in employees_availability:
            for s in range(total_slots):
                model.Add(x[(emp1, s)] + x[(emp2, s)] <= 1)

    # Day-based single shift => generate patterns
    day_pattern_index = {}
    day_worked = {}
    all_patterns = {}

    for e in employees_availability:
        all_patterns[e] = {}
        for d in range(days):
            # Build availability mask for day d
            base = d * 12
            av_mask = []
            for i in range(12):
                slot_id = base + i
                av_mask.append(1 if slot_id in employees_availability[e] else 0)

            valid_pats = generate_day_patterns(av_mask, 6, 8)
            pat_idx = model.NewIntVar(0, len(valid_pats) - 1, f"patIndex_{e}_{d}")
            day_pattern_index[(e, d)] = pat_idx
            day_worked[(e, d)] = model.NewBoolVar(f"dayWorked_{e}_{d}")
            all_patterns[e][d] = valid_pats

    # Enforce only patterns from the table
    for e in employees_availability:
        for d in range(days):
            day_vars = [x[(e, d * 12 + i)] for i in range(12)]
            plus = day_vars + [day_worked[(e, d)], day_pattern_index[(e, d)]]
            table = []
            for idx, pat in enumerate(all_patterns[e][d]):
                row = pat[:12] + [pat[12], idx]
                table.append(row)
            model.AddAllowedAssignments(plus, table)

    # Total shifts => sum(day_worked)
    total_shifts = {}
    diffs = []
    for e in employees_availability:
        sum_days = [day_worked[(e, d)] for d in range(days)]
        total_shifts[e] = model.NewIntVar(0, 7, f"total_shifts_{e}")
        model.Add(sum(sum_days) == total_shifts[e])

    # Desired shifts => if ds in [1,2], forced. else <= ds => minimize difference
    for e in employees_availability:
        ds = desired_shifts[e]
        if ds in [1, 2]:
            model.Add(total_shifts[e] == ds)
        else:
            model.Add(total_shifts[e] <= ds)
            diff_up = model.NewIntVar(0, 7, f"diffUp_{e}")
            diff_down = model.NewIntVar(0, 7, f"diffDown_{e}")
            model.Add(total_shifts[e] - ds <= diff_up)
            model.Add(ds - total_shifts[e] <= diff_down)
            diffs.append(diff_up)
            diffs.append(diff_down)

    # Minimize sum of differences from desired shifts
    model.Minimize(sum(diffs))

    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        # infeasible => None
        return None, None, status, None

    # Build the final schedule
    schedule = {e: [] for e in employees_availability}
    for e in employees_availability:
        for s in range(total_slots):
            if solver.Value(x[(e, s)]) == 1:
                schedule[e].append(s)

    return schedule, solver, status, total_shifts