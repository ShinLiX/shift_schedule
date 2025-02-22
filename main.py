# main.py
"""
Entry point to run the Tkinter GUI for the scheduling application.
"""

import tkinter as tk
from gui_app import SinglePageSchedulerGUI

def main():
    root = tk.Tk()
    app = SinglePageSchedulerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()