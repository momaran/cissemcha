from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from imageio import v2 as imageio  # modern imageio API
from diffusion_mechanism import simulate_pair

# Visual style (close to your previous project)
sns.set_theme(style="darkgrid")
ALPHA = "#e78895"  # α color
BETA  = "#b7b1f2"  # β color


def _cum_from_df(df):
    """
    Given a dataframe with columns:
        col 0 -> 'drained per step'
        col 1 -> 'sourced per step'
    return cumulative arrays: (cum_drained, cum_sourced).
    """
    drained = np.cumsum(df.iloc[:, 0].to_numpy())
    sourced = np.cumsum(df.iloc[:, 1].to_numpy())
    return drained, sourced


def _plot_cumulative(x, y_a, y_b, title, ylabel, ax, color_a=ALPHA, color_b=BETA,
                     xlim=None, ylim=None, legend_title="Spin Type"):
    """Helper to make a tidy cumulative line plot (α/β) like your older style."""
    ax.plot(x, y_a, label="α spins", color=color_a, alpha=0.85, lw=2)
    ax.plot(x, y_b, label="β spins", color=color_b, alpha=0.85, lw=2)
    ax.set_title(title)
    ax.set_xlabel("Time Steps")
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.legend(title=legend_title, loc="upper left")


def run(cfg, outdir: Path,
        q_list=None,
        gif_name: str = "drain_source_qciss.gif",
        frame_dpi: int = 180,
        gif_duration: float = 0.5,
        cleanup_frames: bool = True):
    """
    Produce:
      - drained_vs_steps.png and sourced_vs_steps.png for current cfg.ciss_effect
      - drain_source_qciss.gif sweeping q_CISS

    Parameters
    ----------
    cfg : object
        Your config namespace with fields:
          - n_steps, positions, ciss_effect, voltage_magnitude
          - (optionally) qCISS_vector: iterable of q values for the sweep
    outdir : Path
        Output directory.
    q_list : iterable[float] | None
        If None: use cfg.qCISS_vector if present; else np.linspace(0, 1, 40).
    gif_name : str
        Name of the resulting GIF.
    frame_dpi : int
        DPI for saved frames.
    gif_duration : float
        Seconds per frame in the GIF.
    cleanup_frames : bool
        Delete temporary frames after GIF creation.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------- Static plots for the CURRENT q_CISS ----------
    x = np.arange(cfg.n_steps)
    res0 = simulate_pair(cfg, cfg.voltage_magnitude, trace=False)
    df_a0, df_b0 = res0["df_alpha"], res0["df_beta"]
    aD0, aS0 = _cum_from_df(df_a0)
    bD0, bS0 = _cum_from_df(df_b0)

    y_max_d0 = max(aD0.max(), bD0.max())
    y_max_s0 = max(aS0.max(), bS0.max())

    # Drained vs steps (current q)
    plt.figure(figsize=(6, 4))
    ax = plt.gca()
    _plot_cumulative(
        x, aD0, bD0,
        title=f"Drain Electrons Evolution (q_CISS = {cfg.ciss_effect:.2f})",
        ylabel="Drained Electrons",
        ax=ax,
        xlim=(0, cfg.n_steps),
        ylim=(0, y_max_d0 * 1.05)
    )
    plt.tight_layout()
    plt.savefig(outdir / "drained_vs_steps.png", dpi=300)
    plt.close()

    # Sourced vs steps (current q)
    plt.figure(figsize=(6, 4))
    ax = plt.gca()
    _plot_cumulative(
        x, aS0, bS0,
        title=f"Sourced Electrons (q_CISS = {cfg.ciss_effect:.2f})",
        ylabel="Sourced Electrons",
        ax=ax,
        xlim=(0, cfg.n_steps),
        ylim=(0, y_max_s0 * 1.05)
    )
    plt.tight_layout()
    plt.savefig(outdir / "sourced_vs_steps.png", dpi=300)
    plt.close()

    # ---------- GIF across q_CISS ----------
    if q_list is None:
        q_list = getattr(cfg, "qCISS_vector", None)
    if q_list is None:
        q_list = np.linspace(0.0, 1.0, 40)
    q_list = [float(q) for q in q_list]

    # Cache results per q and compute global y-limits for smooth GIF
    cached = []
    global_ymax_d = 0.0
    global_ymax_s = 0.0

    original_q = float(cfg.ciss_effect)
    try:
        for q in q_list:
            cfg.ciss_effect = q
            res = simulate_pair(cfg, cfg.voltage_magnitude, trace=False)
            df_a, df_b = res["df_alpha"], res["df_beta"]
            aD, aS = _cum_from_df(df_a)
            bD, bS = _cum_from_df(df_b)

            global_ymax_d = max(global_ymax_d, aD.max(), bD.max())
            global_ymax_s = max(global_ymax_s, aS.max(), bS.max())

            cached.append((q, aD, bD, aS, bS))
    finally:
        # Always restore original q
        cfg.ciss_effect = original_q

    # Make frames with fixed limits
    temp_dir = outdir / "_temp_qciss_frames"
    temp_dir.mkdir(parents=True, exist_ok=True)

    frame_paths = []
    for (q, aD, bD, aS, bS) in cached:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
        # Left: Drained cumulative
        _plot_cumulative(
            x, aD, bD,
            title=f"Drained (q_CISS = {q:.2f})",
            ylabel="Cumulative",
            ax=axes[0],
            xlim=(0, cfg.n_steps),
            ylim=(0, global_ymax_d * 1.05),
        )
        # Right: Sourced cumulative
        _plot_cumulative(
            x, aS, bS,
            title=f"Sourced (q_CISS = {q:.2f})",
            ylabel="Cumulative",
            ax=axes[1],
            xlim=(0, cfg.n_steps),
            ylim=(0, global_ymax_s * 1.05),
        )
        fig.tight_layout()
        frame_path = temp_dir / f"frame_qciss_{q:.3f}.png"
        fig.savefig(frame_path, dpi=frame_dpi)
        plt.close(fig)
        frame_paths.append(frame_path)

    # Build GIF
    gif_path = outdir / gif_name
    with imageio.get_writer(gif_path, mode="I", duration=gif_duration) as writer:
        for fp in frame_paths:
            writer.append_data(imageio.imread(fp))

    # Cleanup temporary frames
    if cleanup_frames:
        for fp in frame_paths:
            try:
                os.remove(fp)
            except OSError:
                pass
        try:
            temp_dir.rmdir()
        except OSError:
            # Folder not empty or busy; ignore
            pass

    return {
        "outdir": str(outdir),
        "static_pngs": [
            str(outdir / "drained_vs_steps.png"),
            str(outdir / "sourced_vs_steps.png"),
        ],
        "gif": str(gif_path),
        "n_q_frames": len(frame_paths),
    }