from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from imageio import v2 as imageio  # not used here but handy if you later animate
from diffusion_mechanism import simulate_pair
import matplotlib as mpl

# --- Make sure external LaTeX is NOT used (prevents the crash) ---
mpl.rcParams["text.usetex"] = False
plt.rcParams["font.family"] = "serif"

# Style close to your previous project
sns.set_theme(style="darkgrid")


def _get_title(cfg, current_mode: str) -> str:
    """
    Build a compact, informative title based on cfg fields if available.
    """
    base = f"Current I–V ({current_mode} spins)"
    extras = []
    if hasattr(cfg, "number_spins"):
        extras.append(f"{cfg.number_spins} spins")
    if hasattr(cfg, "temperature"):
        extras.append(f"T = {cfg.temperature:.2f} K")
    if extras:
        return f"{base}, " + ", ".join(extras)
    return base


def _extract_current(res: dict, mode: str) -> float:
    """
    Read current from simulate_pair() result depending on 'mode'.
    Expected keys:
      - 'I_alpha' and/or 'I_beta'; otherwise fallback to 'I_total' if present.
    """
    if mode == "alpha":
        if "I_alpha" in res:
            return float(res["I_alpha"])
        return float(res.get("I_total", 0.0))
    elif mode == "beta":
        if "I_beta" in res:
            return float(res["I_beta"])
        return float(res.get("I_total", 0.0))
    elif mode == "total":
        Ia = float(res.get("I_alpha", 0.0))
        Ib = float(res.get("I_beta", 0.0))
        if (Ia != 0.0) or (Ib != 0.0):
            return Ia + Ib
        return float(res.get("I_total", Ia + Ib))
    else:
        raise ValueError(f"Unknown mode: {mode}")


def run(cfg, outdir: Path,
        mode: str = "alpha",  # 'alpha' | 'beta' | 'total'
        save_csv: bool = False):
    """
    Plot I–V curves for multiple q_CISS values on a single figure.

    Parameters
    ----------
    cfg : object
        Needs: ciss_effect, number_spins, (optional) Temperature, diff_coefficient
        simulate_pair must accept (cfg, voltage, ...) and return a dict with
        'I_alpha'/'I_beta' and/or 'I_total'.
    outdir : Path
        Output directory for the figure (and optional CSV).
    mode : str
        'alpha' (α-only), 'beta' (β-only), or 'total' (sum).
    save_csv : bool
        If True, write a CSV with the I(V) matrix.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Voltages & q sweep
    Vvec = np.array(getattr(cfg, "V_vector",
                            np.linspace(-0.5, 0.5, 101)), dtype=float)
    Qs = list(getattr(cfg, "qCISS_vector",
                      np.linspace(0.0, 1.0, 7)))
    colors = seaborn_color_set3(len(Qs))

    # Prepare storage: rows=V, cols=q
    I_matrix = np.zeros((len(Vvec), len(Qs)), dtype=float)

    original_q = float(cfg.ciss_effect)
    try:
        for i, q in enumerate(Qs):
            cfg.ciss_effect = float(q)

            # Decide spin counts depending on mode
            if mode == "alpha":
                alpha_count = getattr(cfg, "number_spins", 0)
                beta_count = 0
            elif mode == "beta":
                alpha_count = 0
                beta_count = getattr(cfg, "number_spins", 0)
            else:  # 'total'
                # Use both species if your API supports counts;
                # if not, call without counts and rely on defaults.
                alpha_count = getattr(cfg, "number_spins", None)
                beta_count  = getattr(cfg, "number_spins", None)

            for j, V in enumerate(Vvec):
                # Call simulate_pair in a way compatible with your codebase
                if alpha_count is None or beta_count is None:
                    res = simulate_pair(cfg, float(V))
                else:
                    res = simulate_pair(cfg, float(V),
                                        alpha_count=alpha_count,
                                        beta_count=beta_count)
                I_matrix[j, i] = _extract_current(res, mode)
    finally:
        cfg.ciss_effect = original_q  # always restore

    # ----- Plot (ensure LaTeX stays OFF even if globally enabled) -----
    with plt.rc_context({"text.usetex": False}):
        fig, ax = plt.subplots(figsize=(7, 5))
        for i, q in enumerate(Qs):
            # This label uses mathtext (internal) -> no external LaTeX needed
            ax.plot(Vvec, I_matrix[:, i], marker='o', linestyle='-',
                    color=colors[i], label=rf"$q_{{\mathrm{{CISS}}}} = {q:.2f}$")

        title = _get_title(cfg, current_mode=("α-only" if mode == "alpha"
                                              else "β-only" if mode == "beta"
                                              else "total"))
        ax.set_title(title)
        ax.set_xlabel("Voltage (arb. units)")
        ax.set_ylabel("Current (arb. units)")
        ax.grid(True)
        ax.legend()

        # Save figure
        fname = {
            "alpha": "IV_alpha_vs_qciss.png",
            "beta":  "IV_beta_vs_qciss.png",
            "total": "IV_total_vs_qciss.png",
        }[mode]
        fig.savefig(outdir / fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Optional CSV dump
    csv_path = None
    if save_csv:
        import pandas as pd
        df = pd.DataFrame(I_matrix, index=Vvec,
                          columns=[f"qCISS={q:.3f}" for q in Qs])
        df.index.name = "Voltage"
        csv_path = outdir / fname.replace(".png", ".csv")
        df.to_csv(csv_path)

    return {
        "outdir": str(outdir),
        "mode": mode,
        "figure": str(outdir / fname),
        "csv": (str(csv_path) if csv_path else None),
        "q_values": [float(q) for q in Qs],
        "V_values": Vvec.tolist(),
    }


def seaborn_color_set3(n: int):
    """
    Helper to consistently get a Set3 palette with n colors.
    (Avoids re-computing palette in multiple places if you reuse this pattern.)
    """
    return sns.color_palette("Set3", n)
