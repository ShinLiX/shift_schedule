# exporter.py
"""
Module for exporting the schedule to a CSV file.
"""

import csv

def export_schedule_to_csv(
    schedule,
    employees,
    solver,
    total_shifts_vars,
    desired_shifts,
    filename="schedule.csv"
):
    """
    Exports the schedule to a CSV file.

    For each day/hour:
      - occupant1, occupant2 if coverage=2, else occupant2='NA'.
      - If coverage=0 => occupant1= occupant2='NA'.

    Then appends a table of:
      - Employee Name
      - WeeklyHours
      - DesiredShifts
      - ActualShifts
    """
    days = 7
    slots_per_day = 12
    day_names = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]

    # Build (day, hour) -> assigned employees
    day_hour_emps = {d: {h: [] for h in range(slots_per_day)} for d in range(days)}
    for e in employees:
        for slot_id in schedule[e]:
            d = slot_id // slots_per_day
            h = slot_id % slots_per_day
            day_hour_emps[d][h].append(e)

    occupant1 = {(d, h): "NA" for d in range(days) for h in range(slots_per_day)}
    occupant2 = {(d, h): "NA" for d in range(days) for h in range(slots_per_day)}

    # Assign occupant1, occupant2 per day/hour
    for d in range(days):
        for h in range(slots_per_day):
            assigned = sorted(day_hour_emps[d][h])
            if len(assigned) == 1:
                occupant1[(d, h)] = assigned[0]
            elif len(assigned) >= 2:
                occupant1[(d, h)] = assigned[0]
                occupant2[(d, h)] = assigned[1]

    # Create CSV rows => 12 rows, each row is time + 2 columns per day
    headers = ["Time"]
    for d in range(days):
        headers.append(f"{day_names[d]}-1")
        headers.append(f"{day_names[d]}-2")

    rows = []
    for hour in range(slots_per_day):
        time_str = f"{11 + hour}:00 - {12 + hour}:00"
        row_cells = [time_str]
        for d in range(days):
            row_cells.append(occupant1[(d, hour)])
            row_cells.append(occupant2[(d, hour)])
        rows.append(row_cells)

    # Add a blank row
    rows.append([])

    # Add summary table => Name, WeeklyHours, DesiredShifts, ActualShifts
    rows.append(["Name", "WeeklyHours", "DesiredShifts", "ActualShifts"])
    for e in employees:
        total_hours = len(schedule[e])
        wanted = desired_shifts[e]
        got = solver.Value(total_shifts_vars[e])
        rows.append([e, str(total_hours), str(wanted), str(got)])

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row_cells in rows:
            writer.writerow(row_cells)

    print(f"Exported schedule to '{filename}'.")