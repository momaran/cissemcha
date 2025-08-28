# CISS/eMChA Kinetic Monte Carlo (KMC) Simulator

A minimal, modular **kinetic Monte Carlo** code to simulate **spin-selective transport (CISS)** and **magneto-chiral anisotropy (eMChA)** in a 1D channel under voltage bias. The code provides a small GUI to set parameters, several ready-to-run simulations, and figures/animations as outputs.

---

## 0) Quick start

```bash
# Python 3.9–3.12
pip install numpy matplotlib seaborn imageio

# Run (opens the GUI)
python main.py

cissemcha-main/
├── main.py                     # Entry point: opens the GUI and dispatches
├── read_data.py                # Tkinter GUI, validation, derived quantities
├── diffusion_mechanism.py      # Physics + KMC core (Glauber, CISS, eMChA)
├── simulations/
│   ├── Simulation_type.py      # Dispatcher (int → module.run)
│   ├── plot_style.py           # Plot style helpers (matplotlib/seaborn)
│   ├── initial_state.py        # (optional helpers for states)
│   ├── steps_counters/
│   │   └── sim.py              # Drained/Sourced vs steps (+ optional GIF over q_CISS)
│   ├── iv_three_sources/
│   │   └── sim.py              # I–V: unpolarized, alpha-only, beta-only
│   ├── iv_alpha_qciss/
│   │   └── sim.py              # I–V vs multiple q_CISS (α/β/total)
│   └── trajectory_histograms/
│       └── sim.py              # GIF: α and β position histograms along the channel
├── last_config.json            # Auto-saved GUI config (offered on next run)
└── outputs/                    # Results, time-stamped subfolders
