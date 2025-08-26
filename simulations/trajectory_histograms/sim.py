from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import imageio
from diffusion_mechanism import simulate_single

sns.set_theme(style="darkgrid")

ALPHA = "#e78895"
BETA  = "#b7b1f2"


def _counts_at_step(traj, t: int, N: int) -> np.ndarray:
    """
    Return counts per lattice site at time step t.
    Supports several trajectory shapes:
      - ndarray shape (n_spins, n_steps) -> uses column t
      - ndarray shape (n_steps,)         -> single spin over time
      - list/tuple length steps with entry being int or array of positions
    """
    if isinstance(traj, np.ndarray):
        if traj.ndim == 2:
            positions = traj[:, t]
        elif traj.ndim == 1:
            positions = np.array([traj[t]])
        else:
            positions = np.atleast_1d(traj[..., t])
    else:
        positions = np.atleast_1d(traj[t])

    positions = np.asarray(positions, dtype=int)
    # Keep only positions within [0, N] just in case
    positions = positions[(positions >= 0) & (positions <= N)]
    return np.bincount(positions, minlength=N + 1)


def _apply_stem_colors(stem_container, color: str):
    """Set consistent colors for a matplotlib stem plot and hide baseline."""
    markerline, stemlines, baseline = stem_container
    markerline.set_color(color)
    if hasattr(stemlines, "set_color"):
        stemlines.set_color(color)
    else:
        for ln in stemlines:
            ln.set_color(color)
    if baseline is not None:
        baseline.set_visible(False)


def run(cfg, outdir: Path):
    """
    Generate alpha/beta trajectory histograms as animated GIF replicating the old style:
    - stem plots in two stacked subplots
    - per-step counts using np.bincount
    - temporary PNG frames -> single GIF -> cleanup
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Run one simulation
    res = simulate_single(cfg, cfg.voltage_magnitude)
    pa = res["traj_alpha"]
    pb = res["traj_beta"]

    N = int(getattr(cfg, "positions", 0))
    if N <= 0:
        raise ValueError("cfg.positions must be a positive integer")

    def _infer_steps(traj):
        if isinstance(traj, np.ndarray):
            return traj.shape[-1]
        else:
            return len(traj)

    steps = min(_infer_steps(pa), _infer_steps(pb))

    # Old-style temp dir
    temp_dir = outdir / "temp_frames"
    temp_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []

    # Keep every step unless it's huge
    stride = 1 if steps <= 600 else max(1, steps // 600)

    ns = getattr(cfg, "number_spins", None)
    T = getattr(cfg, "temperature", None)
    V = getattr(cfg, "voltage_magnitude", None)
    qc = getattr(cfg, "ciss_effect", None)

    for t in range(0, steps, stride):
        plt.figure(figsize=(10, 5))

        # α subplot
        plt.subplot(2, 1, 1)
        alpha_counts = _counts_at_step(pa, t, N)
        stem_a = plt.stem(np.arange(len(alpha_counts)),
                          alpha_counts,
                          linefmt='-',
                          markerfmt='o',
                          basefmt=" ")
        _apply_stem_colors(stem_a, ALPHA)
        plt.xlabel("Position")
        plt.ylabel("Number of α spins")
        if ns is not None and T is not None and V is not None and qc is not None:
            plt.title(
                f"Distribution of α spins. Nstep {t} for {ns} total spins, "
                f"D = {T:.2f}, V = {V:.1f}, qciss = {qc:.2f}"
            )
        else:
            plt.title(f"Distribution of α spins. Nstep {t}")

        # β subplot
        plt.subplot(2, 1, 2)
        beta_counts = _counts_at_step(pb, t, N)
        stem_b = plt.stem(np.arange(len(beta_counts)),
                          beta_counts,
                          linefmt='-',
                          markerfmt='o',
                          basefmt=" ")
        _apply_stem_colors(stem_b, BETA)
        plt.xlabel("Position")
        plt.ylabel("Number of β spins")
        if ns is not None and T is not None and V is not None and qc is not None:
            plt.title(
                f"Distribution of β spins. Nstep {t} for {ns} total spins, "
                f"T = {T:.2f}, V = {V:.1f}, qciss = {qc:.2f}"
            )
        else:
            plt.title(f"Distribution of β spins. Nstep {t}")

        plt.tight_layout()
        frame_path = temp_dir / f"frame_{t:05d}.png"
        plt.savefig(frame_path, dpi=200)
        plt.close()
        frame_paths.append(str(frame_path))

    gif_path = outdir / "spins_evolution.gif"
    with imageio.get_writer(gif_path, mode='I', duration=0.1) as writer:
        for fp in frame_paths:
            writer.append_data(imageio.imread(fp))

    for fp in frame_paths:
        try:
            os.remove(fp)
        except OSError:
            pass
    try:
        os.rmdir(temp_dir)
    except OSError:
        pass

    return {
        "outdir": str(outdir),
        "gif": str(gif_path),
        "frames_used": len(frame_paths),
        "stride": stride,
        "steps": steps
    }

