import tkinter as tk
from tkinter import filedialog, messagebox

from .simulation import run_simulation


def start_gui():
    """Launch a minimal GUI to run the diffusion simulation."""
    root = tk.Tk()
    root.title("CISS-eMChA Simulation")

    text = tk.Text(root, width=80, height=20)
    text.pack(padx=10, pady=10)

    def run():
        path = filedialog.askopenfilename(
            title="Select configuration file",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            summary, _, _ = run_simulation(path)
        except Exception as exc:
            messagebox.showerror("Error", str(exc))
            return

        text.delete("1.0", tk.END)
        text.insert(tk.END, summary.to_string(index=False))

    tk.Button(root, text="Run Simulation", command=run).pack(pady=5)
    root.mainloop()
