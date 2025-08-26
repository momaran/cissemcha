# read_data.py
import tkinter as tk
from tkinter import ttk, messagebox
import json
import os
import numpy as np

# Import the SimulationType enum from our dispatcher
from simulations.Simulation_type import SimulationType

CONFIG_FILE = os.path.join(os.path.dirname(__file__), "last_config.json")

def _cast_value(raw: str, caster):
    raw = raw.strip()
    if caster is int:
        return int(float(raw))  # allow '10.0' in an int field
    return caster(raw)

def _load_previous():
    if not os.path.exists(CONFIG_FILE):
        return None
    try:
        with open(CONFIG_FILE, "r") as f:
            data = json.load(f)
        # best-effort cast back to numbers where it makes sense
        for k, v in list(data.items()):
            if isinstance(v, str):
                try:
                    if "." in v or "e" in v.lower():
                        data[k] = float(v)
                    else:
                        data[k] = int(v)
                except Exception:
                    pass
        return data
    except Exception:
        return None

def _save_config(cfg: dict):
    with open(CONFIG_FILE, "w") as f:
        json.dump(cfg, f, indent=2)

def read_data():
    # Ask whether to reuse last config if present
    prev = _load_previous()
    if prev is not None:
        reuse = messagebox.askyesno("Reuse Previous Values?",
                                    "Do you want to reuse the previous values?")
        if reuse:
            return prev

    root = tk.Tk()
    root.title("CISS / eMChA Monte Carlo • Parameters")
    root.geometry("680x820")
    root.configure(bg="#DDEEFF")

    entries = {}

    # --- Simulation selector
    ttk.Label(root, text="Simulation Type:").pack(pady=(10, 0))
    sim_type_var = tk.StringVar()

    sim_choices = [e.name for e in SimulationType]
    sim_combo = ttk.Combobox(root, values=sim_choices, textvariable=sim_type_var, state="readonly")
    sim_combo.current(0)
    sim_combo.pack()

    # --- Numeric parameters (default values chosen to be sensible)
    params = [
        ("Temperature (K):",            "Temperature",           "300.0", float),
        ("Activation Energy U_eff (eV):","activation_energy",     "0.05",  float),
        ("Attempt Time τ0 (s):",        "tau_0",                 "1e-12", float),
        ("Molecule Length L (m):",      "molecule_length",       "5e-9",  float),
        ("Positions (sites):",          "positions",             "15",    int),
        ("Step size dx (m):",           "dx",                    "1e-10", float),
        ("Number of steps:",            "n_steps",               "400",   int),
        ("Number of electrons:",        "number_spins",          "400",   int),
        ("Spin ratio α/β:",             "spin_ratio",            "1",     int),
        ("Voltage amplitude (V):",      "voltage_magnitude",     "0.1",   float),
        ("Voltage frequency (Hz):",     "voltage_freq",          "0.0",   float),
        ("CISS coupling q_ciss [0–1]:", "ciss_effect",           "0.5",   float),
        ("eMChA coupling β_emcha:",     "beta_emcha",            "0.0",   float),
        ("External B field (T):",       "B_ext",                 "0.0",   float),
        ("Chirality (+1 / -1):",        "chirality",             "1",     int),
        ("Jump probability p_jump:",    "jump_probability",      "0.15",  float),
        ("GIF frames (if used):",       "gif_frames",            "100",   int),
    ]

    # Grid container
    grid = ttk.Frame(root)
    grid.pack(fill="x", pady=10, padx=14)

    for i, (label, key, default, caster) in enumerate(params):
        fr = ttk.Frame(grid)
        fr.grid(row=i, column=0, sticky="ew", pady=3)
        ttk.Label(fr, text=label, width=32).pack(side="left")
        ent = ttk.Entry(fr)
        ent.insert(0, default)
        ent.pack(side="left", fill="x", expand=True)
        entries[key] = (ent, caster)

    # Toggles
    toggles = ttk.Frame(root)
    toggles.pack(fill="x", pady=8, padx=14)
    use_ciss_var = tk.BooleanVar(value=True)
    ttk.Checkbutton(toggles, text="Enable CISS", variable=use_ciss_var).pack(side="left", padx=4)
    use_emcha_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(toggles, text="Enable eMChA", variable=use_emcha_var).pack(side="left", padx=4)
    self_consistent_B_var = tk.BooleanVar(value=False)
    ttk.Checkbutton(toggles, text="Self-consistent internal B(I)", variable=self_consistent_B_var).pack(side="left", padx=4)

    # Optional helix pitch for self-consistent internal field μ0 I / pitch
    pitch_fr = ttk.Frame(root)
    pitch_fr.pack(fill="x", pady=2, padx=14)
    ttk.Label(pitch_fr, text="Helix pitch (m) for B_int = μ0 I / pitch: ", width=32).pack(side="left")
    pitch_entry = ttk.Entry(pitch_fr)
    pitch_entry.insert(0, "0.0")
    pitch_entry.pack(side="left", fill="x", expand=True)

    def on_help():
        help_text = (
            "SIMULATIONS\n"
            "1) STEPS_COUNTERS: Drain/Sourced vs steps; optional GIF across q_ciss.\n"
            "2) IV_THREE_SOURCES: I–V for unpolarized, beta-only, alpha-only.\n"
            "3) IV_ALPHA_QCISS: I–V for alpha electrons sweeping q_ciss.\n"
            "4) TRAJECTORY_HISTOGRAMS: GIF of α and β electron position histograms along the channel.\n\n"
            "NOTES\n"
            "- dt is derived from τ = τ0 exp(U_eff/kBT) and p_jump: dt = p_jump / Γ, Γ≈2/τ.\n"
            "- For self-consistent internal field, set Helix pitch > 0 and enable the toggle.\n"
        )
        top = tk.Toplevel(root)
        top.title("Help")
        txt = tk.Text(top, wrap="word", height=22, width=72)
        txt.insert("1.0", help_text)
        txt.config(state="disabled")
        txt.pack(padx=10, pady=10)

    def on_ok():
        try:
            cfg = {}
            # sim type
            name = sim_type_var.get()
            enum_val = SimulationType[name].value
            cfg["simulation_type"] = enum_val

            # numerics
            for key, (ent, caster) in entries.items():
                cfg[key] = _cast_value(ent.get(), caster)

            # toggles
            cfg["use_ciss"] = bool(use_ciss_var.get())
            cfg["use_emcha"] = bool(use_emcha_var.get())
            cfg["self_consistent_B"] = bool(self_consistent_B_var.get())

            # helix pitch
            try:
                cfg["helix_pitch"] = float(pitch_entry.get().strip())
            except Exception:
                cfg["helix_pitch"] = 0.0

            # Derived constants (SI + eV)
            kB_eV = 8.617333262145e-5  # eV/K
            q = 1.602176634e-19         # C
            cfg["Boltzmann_constant"] = kB_eV
            cfg["electron_charge"] = q

            # Attempt time and rate, dt consistent with p_jump (kinetic MC logic)
            tau = cfg["tau_0"] * np.exp(cfg["activation_energy"] / (cfg["Temperature"] * kB_eV))
            Gamma_pair = 2.0 / tau  # total attempt rate (left+right)
            dt = max(1e-18, cfg["jump_probability"]) / Gamma_pair
            cfg["relaxation_time"] = tau
            cfg["relaxation_rate"] = 1.0 / tau
            cfg["dt"] = dt

            # Diffusion coefficient (Einstein-like estimate for 1D hop)
            cfg["diff_coefficient"] = (cfg["dx"] ** 2) / (2.0 * tau)

            # Field per site for convenience
            cfg["dV"] = cfg["voltage_magnitude"] / max(cfg["positions"], 1)

            # Spin populations
            ns = int(cfg["number_spins"])
            sr = max(0, int(cfg["spin_ratio"]))
            alpha = ns // (sr + 1)
            beta = ns - alpha
            cfg["alpha_spins"] = int(alpha)
            cfg["beta_spins"] = int(beta)

            _save_config(cfg)
            root.destroy()
            # return via closure trick: attach to function attribute
            read_data._result = cfg
        except Exception as e:
            messagebox.showerror("Error", str(e))

    btns = ttk.Frame(root); btns.pack(pady=10)
    ttk.Button(btns, text="Help", command=on_help).pack(side="left", padx=6)
    ttk.Button(btns, text="OK", command=on_ok).pack(side="left", padx=6)

    root.mainloop()
    return getattr(read_data, "_result", prev or {})
