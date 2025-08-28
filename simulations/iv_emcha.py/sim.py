# iv_emcha/sim.py
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class EMCHAConfig:
    r0: float = 1.0           # Base resistance at B=0 [Ohm]
    alpha_B2: float = 0.0     # Even-in-B coefficient (set >0 to break tangency)
    beta_BI: float = 0.25     # Odd-in-(B*I) coefficient (controls nonreciprocity)
    chi: int = +1             # Molecular chirality (+1 right-handed, -1 left-handed)
    spin_channel: str = "alpha"  # For labeling only
    v_min: float = -0.2
    v_max: float = 0.2
    n_points: int = 401
    b_values: tuple[float, float, float] = (+1.0, 0.0, -1.0)  # relative B (arb. units)
    include_B2: bool = False  # if True, add alpha_B2 * B^2 term (no longer tangent)
    outdir: Path = Path("outputs")


def current_from_voltage_emcha(V: np.ndarray, B: float, cfg: EMCHAConfig) -> np.ndarray:
    """
    Phenomenological EMCA model with:
        R(I, B) = R0 * ( 1 + alpha*B^2 + beta * chi * B * I )
    => V = I * R(I, B) = R0*(1 + alpha*B^2)*I + R0*(beta*chi*B)*I^2

    For each V, solve the quadratic in I:
        C I^2 + A I - V = 0
      where A = R0*(1 + alpha*B^2), C = R0*beta*chi*B

    Select the root that matches the small-V expansion I ~ V/A.
    """
    A = cfg.r0 * (1.0 + (cfg.alpha_B2 * B * B if cfg.include_B2 else 0.0))
    C = cfg.r0 * cfg.beta_BI * cfg.chi * B

    # Handle C ≈ 0 (B=0 or beta=0) → strictly Ohmic
    if np.isclose(C, 0.0):
        return V / A

    # Quadratic solution: pick the root continuous with I ~ V/A at small V
    # I = (-A + s*sqrt(A^2 + 4*C*V)) / (2*C), with s = sign(V) for continuity
    D = A * A + 4.0 * C * V
    # Ensure numerical safety for tiny negative D from float errors
    D = np.where(D < 0, 0.0, D)
    s = np.sign(V + 1e-15)
    I = (-A + s * np.sqrt(D)) / (2.0 * C)
    return I


def simulate(cfg: EMCHAConfig):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = cfg.outdir / f"iv_emcha_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)

    V = np.linspace(cfg.v_min, cfg.v_max, cfg.n_points)
    curves = []
    for B in cfg.b_values:
        I = current_from_voltage_emcha(V, B, cfg)
        label = f"B={B:+g} ({cfg.spin_channel})"
        curves.append((B, label, V, I))
        df = pd.DataFrame({"V": V, "I": I})
        df.to_csv(run_dir / f"IV_emcha_B_{B:+g}.csv", index=False)

    # Plot
    plt.figure(figsize=(6.0, 4.2), dpi=140)
    for _, label, V, I in curves:
        plt.plot(V, I, lw=2, label=label)
    title = "I–V with EMCA (odd BI term{})".format(
        " + even B² term" if cfg.include_B2 else ""
    )
    plt.title(title)
    plt.xlabel("Voltage V (arb.)")
    plt.ylabel(f"Current I (alpha channel)")
    plt.legend()
    plt.grid(True, alpha=0.35)
    png_path = run_dir / "IV_emcha.png"
    plt.tight_layout()
    plt.savefig(png_path)
    print(f"[OK] Saved figure: {png_path}")
    print(f"[OK] Saved CSVs in: {run_dir}")


def main():
    p = argparse.ArgumentParser(
        description="Simulate I–V curves with electrical magnetochiral anisotropy (EMCA)."
    )
    p.add_argument("--r0", type=float, default=1.0, help="Base resistance at B=0 [Ohm]")
    p.add_argument("--beta_BI", type=float, default=0.25, help="Odd BI coefficient")
    p.add_argument("--alpha_B2", type=float, default=0.05, help="Even B^2 coefficient")
    p.add_argument("--include_B2", action="store_true", help="Include even-in-B term")
    p.add_argument("--chi", type=int, default=+1, choices=[-1, +1], help="Chirality")
    p.add_argument("--spin", type=str, default="alpha", help="Spin channel label")
    p.add_argument("--vmin", type=float, default=-0.2)
    p.add_argument("--vmax", type=float, default=0.2)
    p.add_argument("--points", type=int, default=401)
    p.add_argument("--Bvals", type=float, nargs="+", default=[+1.0, 0.0, -1.0],
                   help="List of B values (arb. units)")
    p.add_argument("--out", type=str, default="outputs")
    args = p.parse_args()

    cfg = EMCHAConfig(
        r0=args.r0,
        alpha_B2=args.alpha_B2,
        beta_BI=args.beta_BI,
        chi=args.chi,
        spin_channel=args.spin,
        v_min=args.vmin,
        v_max=args.vmax,
        n_points=args.points,
        b_values=tuple(args.Bvals),
        include_B2=args.include_B2,
        outdir=Path(args.out),
    )
    simulate(cfg)


if __name__ == "__main__":
    main()
