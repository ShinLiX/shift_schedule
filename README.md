# **Employee Shift Scheduling Tool**

## **Overview**
This Python-based scheduling tool assigns employees to shifts while satisfying various hard and soft constraints, ensuring optimal shift coverage and fairness in employee allocation. The tool is powered by **Google OR-Tools (CP-SAT Solver)** and is designed to handle complex scheduling requirements efficiently. 

The generated schedule is exported to a CSV with intuitive tables for easy analysis and usage.

---

## **Features**

### **1. Constraint Satisfaction**
The tool enforces the following **hard constraints**:
- **Coverage Requirement**: Each time slot must have exactly 2 employees.
- **Shift Lengths**: Employees are assigned single shifts of 6–8 hours (or shorter if daily availability is less than 6 hours).
- **Weekly Hour Limits**: Employees cannot exceed their predefined weekly working hours.
- **Forbidden Pairs**: Prevents specific pairs of employees from being assigned to the same shift.
- **Shift Continuity**: Day-based single-shift patterns are strictly enforced.

### **2. Soft Objectives**
The tool minimizes deviations from employees' preferences, including:
- Assigning employees their desired number of shifts whenever possible.
- Ensuring employees requesting exactly 1 or 2 shifts are strictly matched.

### **3. CSV Export**
The resulting schedule is exported in a user-friendly format:
1. **Time vs. Day Grid**: Displays shifts in a grid format, preserving occupant continuity by listing each employee’s name at the top of their contiguous block.
2. **Summary Table**: Includes details for each employee:
   - Weekly hours
   - Desired shift count
   - Actual shift count

---

## **Technologies Used**
- **Python**
- **Google OR-Tools (CP-SAT Solver)** for constraint satisfaction modeling
- **pandas** for data manipulation and CSV export

---

## **How It Works**
1. Define employees, availability, and constraints in an input configuration file (e.g., JSON or CSV).
2. The tool constructs a **constraint satisfaction model** with the following logic:
   - Coverage requirements per slot (`AddEquality`).
   - Single-shift patterns and continuity (`AddAllowedAssignments`).
   - Weekly limits and forbidden pairs.
3. The model is solved using the **CP-SAT solver**, optimizing for fairness and minimizing deviations from employee preferences.
4. Export the schedule to a CSV file containing:
   - A **time vs. day grid** for shift visualization.
   - A **summary table** for performance review and adjustments.

---

## **Getting Started**

### **Prerequisites**
- Python 3.8+
- Google OR-Tools
- pandas

### **Installation**
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/shift-scheduler.git
   cd shift-scheduler