import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess

def run_script():
    database_file = database_entry.get()
    evidence_file = evidence_entry.get()
    suspect_file = suspect_entry.get()
    threads = threads_entry.get()
    log10_flag = log10_var.get()

    if not database_file or not evidence_file or not suspect_file or not threads:
        messagebox.showerror("Input Error", "All fields are required.")
        return

    command = ["python", "pRMNE.py", "-d", database_file, "-e", evidence_file, "-s", suspect_file, "-t", threads]
    if log10_flag:
        command.append("-l")

    try:
        result = subprocess.run(command, capture_output=True, text=True)
        output_text.delete(1.0, tk.END)
        output_text.insert(tk.END, result.stdout)
    except Exception as e:
        messagebox.showerror("Execution Error", str(e))

def browse_file(entry):
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)

# Create the main window
root = tk.Tk()
root.title("pRMNE GUI")

# Configure grid to expand with window size
root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(7, weight=1)

# Add software name label
tk.Label(root, text="The probability of Random Man Not Excluding (pRMNE) V1.0.0").grid(row=0, column=0, columnspan=3, padx=10, pady=10)

# Create and place the labels and entries
tk.Label(root, text="Database File:").grid(row=1, column=0, padx=10, pady=5, sticky=tk.W)
database_entry = tk.Entry(root, width=50)
database_entry.grid(row=1, column=1, padx=10, pady=5, sticky=tk.EW)
tk.Button(root, text="Browse", command=lambda: browse_file(database_entry)).grid(row=1, column=2, padx=10, pady=5)

tk.Label(root, text="Evidence File:").grid(row=2, column=0, padx=10, pady=5, sticky=tk.W)
evidence_entry = tk.Entry(root, width=50)
evidence_entry.grid(row=2, column=1, padx=10, pady=5, sticky=tk.EW)
tk.Button(root, text="Browse", command=lambda: browse_file(evidence_entry)).grid(row=2, column=2, padx=10, pady=5)

tk.Label(root, text="Suspect File:").grid(row=3, column=0, padx=10, pady=5, sticky=tk.W)
suspect_entry = tk.Entry(root, width=50)
suspect_entry.grid(row=3, column=1, padx=10, pady=5, sticky=tk.EW)
tk.Button(root, text="Browse", command=lambda: browse_file(suspect_entry)).grid(row=3, column=2, padx=10, pady=5)

tk.Label(root, text="Threads:").grid(row=4, column=0, padx=10, pady=5, sticky=tk.W)
threads_entry = tk.Entry(root, width=50)
threads_entry.insert(0, "2")  # Set default value to 2
threads_entry.grid(row=4, column=1, padx=10, pady=5, sticky=tk.EW)

log10_var = tk.BooleanVar()
tk.Checkbutton(root, text="Log10 Conversion", variable=log10_var).grid(row=5, column=1, padx=10, pady=5, sticky=tk.W)

tk.Button(root, text="Run", command=run_script).grid(row=6, column=1, padx=10, pady=20)

output_text = tk.Text(root, height=15, width=80)
output_text.grid(row=7, column=0, columnspan=3, padx=10, pady=5, sticky=tk.NSEW)

# Start the main event loop
root.mainloop()