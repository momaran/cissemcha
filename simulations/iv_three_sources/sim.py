from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from diffusion_mechanism import simulate_pair

sns.set_theme(style="darkgrid")
ALPHA = "#e78895"
BETA  = "#b7b1f2"

def run(cfg, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    Vvec = np.linspace(-10.0, 10.0, 100, dtype=float)

    I_unpol = []; I_alpha = []; I_beta = []

    for V in Vvec:
        res = simulate_pair(cfg, float(V))
        I_unpol.append(res["I_alpha"] + res["I_beta"])

    for V in Vvec:
        res = simulate_pair(cfg, float(V), alpha_count=cfg.number_spins, beta_count=0)
        I_alpha.append(res["I_alpha"])

    for V in Vvec:
        res = simulate_pair(cfg, float(V), alpha_count=0, beta_count=cfg.number_spins)
        I_beta.append(res["I_beta"])

    plt.figure(figsize=(10,6))
    plt.plot(Vvec, I_unpol, label="Unpolarized", color="black")
    plt.plot(Vvec, I_alpha, label="α only", color=ALPHA)
    plt.plot(Vvec, I_beta,  label="β only", color=BETA)
    title = (f"I–V | T={cfg.Temperature:.0f}K, q_CISS={cfg.ciss_effect:.2f}, "
             f"beta_eMChA={getattr(cfg,'beta_emcha',0.0):.2e}, "
             f"Ueff={cfg.activation_energy:.2f} eV, tau0={cfg.tau_0:.1e}s")
    plt.title(title)
    plt.xlabel("Voltage (V)"); plt.ylabel("Current (A)")
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir/"IV_three_cases.png", dpi=300); plt.close()
    return {"outdir": str(outdir)}
