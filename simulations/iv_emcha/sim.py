# simulations/iv_emcha/sim.py
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from diffusion_mechanism import simulate_pair

sns.set_theme(style="darkgrid")
ALPHA = "#e78895"

def _clone_cfg(cfg, **updates):
    """Make a shallow copy of cfg (SimpleNamespace) with updated fields."""
    d = dict(vars(cfg))
    d.update(updates)
    return SimpleNamespace(**d)

def run(cfg, outdir: Path):
    """
    I–V for alpha channel with eMChA:
    - Three curves: B=+|B_ext|, B=0, B=-|B_ext|.
    - Forces eMChA ON and CISS OFF to keep tangency at V=0 (no even-in-B term).
    - Uses only keys already provided by read_data.py.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Voltage sweep centered at 0 using your GUI amplitude
    Vmax = float(getattr(cfg, "voltage_magnitude", 0.1))
    NPTS = 201
    Vvec = np.linspace(-Vmax, +Vmax, NPTS, dtype=float)

    # B values: use |B_ext| from GUI; if 0, fall back to 0.5 T to visualize the effect
    B0_gui = float(getattr(cfg, "B_ext", 0.0))
    B0 = abs(B0_gui) if abs(B0_gui) > 0.0 else 0.5
    B_list = [ +B0, 0.0, -B0 ]

    # Spin counts: alpha-only as you asked
    n_alpha = int(getattr(cfg, "number_spins", 400))
    n_beta  = 0

    # eMChA strength from GUI
    beta_emcha = float(getattr(cfg, "beta_emcha", 0.0))

    # Collect curves
    curves = []
    for B in B_list:
        I_alpha = []
        # Force the physics you want for this sim:
        # - eMChA ON, CISS OFF
        # - set external B for this branch
        cfgB = _clone_cfg(
            cfg,
            B_ext=B,
            use_emcha=True,
            enable_emcha=True,
            use_ciss=False,
            enable_ciss=False,
            beta_emcha=beta_emcha,
        )
        for V in Vvec:
            res = simulate_pair(cfgB, float(V), alpha_count=n_alpha, beta_count=n_beta, trace=False)
            I_alpha.append(res["I_alpha"])
        curves.append((B, np.array(I_alpha, dtype=float)))

    # ---- Plot
    plt.figure(figsize=(10, 6))
    plt.plot(Vvec, curves[0][1], label=f"B=+{B0:g} T (α)", color=ALPHA, lw=2)
    plt.plot(Vvec, curves[1][1], label="B=0 (α)", color="black", lw=2, ls="--")
    plt.plot(Vvec, curves[2][1], label=f"B=-{B0:g} T (α)", color="#4444aa", lw=2)

    title = (f"I–V (alpha) • eMChA only | "
             f"T={getattr(cfg,'Temperature',300):.0f}K, "
             f"β_emcha={beta_emcha:.2e} eV/(T·A), "
             f"U_eff={getattr(cfg,'activation_energy',0.05):.2f} eV, "
             f"τ0={getattr(cfg,'tau_0',1e-12):.1e}s")
    plt.title(title)
    plt.xlabel("Voltage (V)"); plt.ylabel("Current I_alpha (A)")
    plt.legend(); plt.tight_layout()
    fig_path = outdir / "IV_emcha_alpha.png"
    plt.savefig(fig_path, dpi=300); plt.close()

    # ---- CSV
    import pandas as pd
    df = pd.DataFrame({
        "V": Vvec,
        f"I_alpha_B+{B0:g}": curves[0][1],
        "I_alpha_B0":        curves[1][1],
        f"I_alpha_B-{B0:g}": curves[2][1],
    })
    csv_path = outdir / "IV_emcha_alpha.csv"
    df.to_csv(csv_path, index=False)

    print(f"[OK] Saved: {fig_path}")
    print(f"[OK] Saved: {csv_path}")
    return {"outdir": str(outdir)}
