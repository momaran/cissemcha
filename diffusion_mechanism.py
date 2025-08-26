from __future__ import annotations
import numpy as np
import pandas as pd

E_CHARGE = 1.602176634e-19  # C
K_B_J    = 1.380649e-23     # J/K
K_B_eV   = 8.617333262145e-5  # eV/K
MU0      = 4e-7 * np.pi       # H/m

def attempt_pair_rate(T, U_eff_eV, tau0):
    """Arrhenius attempt time -> total (left+right) rate Γ_pair ~ 2/τ."""
    tau = tau0 * np.exp(U_eff_eV / (K_B_eV * T))
    Gamma_pair = 2.0 / tau
    return Gamma_pair, tau

def field_energy_bias_eV(V, L, positions, T):
    """
    Bias per hop (right-left) from uniform E-field over one lattice spacing a = L/(N-1).
    Positive ΔE favors right hops in Glauber splitting.
    """
    a = L / max(positions - 1, 1)
    E_field = V / max(L, 1e-12)
    dE_J = E_CHARGE * E_field * a
    return dE_J / E_CHARGE  # eV

def glauber_rates(DeltaE_eV, T, Gamma_pair):
    """
    Glauber splitting of total rate into right/left; implies tanh-like asymmetry:
      wR - wL = Gamma_pair * tanh(DeltaE/(2kT)).
    """
    x = DeltaE_eV / (K_B_eV * T)
    if x > 50:
        return Gamma_pair, 0.0
    if x < -50:
        return 0.0, Gamma_pair
    wR = Gamma_pair / (1.0 + np.exp(-x))
    wL = Gamma_pair - wR
    return float(wR), float(wL)

# --- New: CISS as a target spin polarization P(V) enforced via artanh ---
def ciss_energy_split_eV(T, V, chirality, P0=None, Pmax=None, q=None):
    """
    Return spin-independent CISS energy scale ΔE_ciss (apply ± with spin later).
    We enforce: tanh(ΔE_ciss/(2kT)) = chirality * P(V).
    Options:
      - Constant polarization: P(V) = P0 in (0,1)
      - Saturating vs bias:   P(V) = Pmax * tanh(q * V), with 0 < Pmax < 1
    """
    kT_eV = K_B_eV * T
    if P0 is not None:
        P = float(P0)
    else:
        if Pmax is None or q is None:
            return 0.0
        P = float(Pmax) * np.tanh(float(q) * V)
    P = np.clip(chirality * P, -0.999, 0.999)
    return 2.0 * kT_eV * np.arctanh(P)  # eV

# --- New: eMChA as a chirality-odd, spin-independent energetic bias ---
def emcha_bias_eV(B_parallel, I_est, beta_emcha, chirality):
    """
    Linear eMChA energetic bias: ΔE_em = chirality * beta_emcha * (B · I).
    Units: beta_emcha in eV / (T·A); B in T; I in A; result in eV.
    Odd in chirality, B, and I. Spin-independent.
    """
    if beta_emcha == 0.0 or B_parallel == 0.0 or I_est == 0.0:
        return 0.0
    return chirality * float(beta_emcha) * (B_parallel * I_est)

def run_channel(cfg, V, Gamma_pair, dt, spin_sigma, Np,
                trace_positions=False, single_particle=False):
    """
    One spin channel.
    * spin_sigma: +1 for α, -1 for β (consistent with labeling below)
    * Periodic boundary conditions (circular):
      - Right jump from N-1 -> counts DRAIN and reappears at 0
      - Left  jump from 0   -> counts SOURCE and reappears at N-1
    """
    N = int(cfg.positions)
    if single_particle:
        Np = 1
    pos = np.random.randint(0, N, size=Np)

    drained = np.zeros(cfg.n_steps, dtype=int)
    sourced = np.zeros(cfg.n_steps, dtype=int)
    hist = np.zeros((cfg.n_steps, N), dtype=int) if trace_positions and not single_particle else None
    traj = np.zeros(cfg.n_steps, dtype=int) if single_particle else None

    enable_ciss  = bool(getattr(cfg, "enable_ciss",  getattr(cfg, "use_ciss", True)))
    enable_emcha = bool(getattr(cfg, "enable_emcha", getattr(cfg, "use_emcha", False)))

    # Field bias per hop (independent of spin)
    DeltaE_field = field_energy_bias_eV(V, cfg.molecule_length, N, cfg.Temperature)

    # CISS energy base (independent of spin; sign applied with spin_sigma)
    DeltaE_ciss_base = 0.0
    if enable_ciss:
        P0_cfg   = getattr(cfg, "ciss_P0", None)
        Pmax_cfg = getattr(cfg, "ciss_Pmax", None)
        q_cfg    = getattr(cfg, "ciss_q", None)
        # Backward-compat: if nothing provided, use ciss_effect as constant P0
        if P0_cfg is None and Pmax_cfg is None:
            P0_cfg = float(getattr(cfg, "ciss_effect", 0.0))
        DeltaE_ciss_base = ciss_energy_split_eV(
            cfg.Temperature, V, getattr(cfg, "chirality", 1),
            P0=P0_cfg, Pmax=Pmax_cfg, q=q_cfg
        )

    I_est = 0.0  # cumulative time-averaged current estimate (A) for this channel
    for t in range(cfg.n_steps):
        if hist is not None:
            hist[t] = np.bincount(pos, minlength=N)
        if traj is not None:
            traj[t] = pos[0]

        # Total B field (optional self-consistent contribution from I_est)
        B_tot = float(getattr(cfg, "B_ext", 0.0))
        if getattr(cfg, "self_consistent_B", False) and float(getattr(cfg, "helix_pitch", 0.0)) > 0.0:
            B_tot += MU0 * I_est / float(cfg.helix_pitch)

        # Spin-independent eMChA bias (uses current estimate up to previous steps)
        beta_emcha = float(getattr(cfg, "beta_emcha", 0.0)) if enable_emcha else 0.0
        DeltaE_em = emcha_bias_eV(B_tot, I_est, beta_emcha, getattr(cfg, "chirality", 1))

        # Total energy bias for this spin
        DeltaE_total = DeltaE_field + (spin_sigma * DeltaE_ciss_base) + DeltaE_em

        # KMC step (Glauber)
        wR, wL = glauber_rates(DeltaE_total, cfg.Temperature, Gamma_pair)
        pR = min(0.9999, wR * dt)
        pL = min(0.9999, wL * dt)

        prev = pos.copy()
        r = np.random.rand(Np)
        moveR = (r < pR)
        r2 = np.random.rand(Np)
        moveL = (~moveR) & (r2 < pL)

        # --- periodic boundary updates using prev state ---
        right_edge = moveR & (prev == N - 1)  # crosses DRAIN, wrap to 0
        left_edge  = moveL & (prev == 0)      # crosses SOURCE, wrap to N-1
        insideR    = moveR & (prev < N - 1)
        insideL    = moveL & (prev > 0)

        pos = prev  # start from prev and update
        pos[insideR] += 1
        pos[insideL] -= 1
        pos[right_edge] = 0
        pos[left_edge]  = N - 1

        drained_now = int(np.sum(right_edge))
        sourced_now = int(np.sum(left_edge))

        drained[t] = drained_now
        sourced[t] = sourced_now

        # Update current estimate for next step (cumulative average)
        Q = E_CHARGE * np.sum(drained[:t + 1] - sourced[:t + 1])
        T_sofar = (t + 1) * dt
        I_est = Q / max(T_sofar, 1e-15)

    label = 'α' if spin_sigma > 0 else 'β'
    df = pd.DataFrame({
        f"Drained {label} spins": drained,
        f"Sourced {label} spins": sourced,
    })
    return df, hist, traj

def simulate_pair(cfg, V, alpha_count=None, beta_count=None, trace=False):
    """Simulate both spin channels independently and return currents & traces."""
    Gamma_pair, tau = attempt_pair_rate(
        cfg.Temperature,
        float(getattr(cfg, "activation_energy", 0.05)),
        cfg.tau_0
    )

    # jump_probability preferred; fallback to diff_coefficient for backward-compat
    p_jump = float(getattr(cfg, "jump_probability", getattr(cfg, "diff_coefficient", 0.15)))
    p_jump = max(1e-4, min(0.3, p_jump))
    dt = p_jump / Gamma_pair

    n_alpha = int(getattr(cfg, "alpha_spins", getattr(cfg, "number_spins", 1000) // 2)) if alpha_count is None else int(alpha_count)
    n_beta  = int(getattr(cfg, "beta_spins",  getattr(cfg, "number_spins", 1000) - n_alpha)) if beta_count  is None else int(beta_count)

    df_alpha, hist_a, traj_a = run_channel(cfg, V, Gamma_pair, dt, spin_sigma=+1, Np=n_alpha, trace_positions=trace)
    df_beta,  hist_b, traj_b = run_channel(cfg, V, Gamma_pair, dt, spin_sigma=-1, Np=n_beta,  trace_positions=trace)

    # Net charge that crossed drain minus source in total time
    Q_alpha = E_CHARGE * (df_alpha.iloc[:, 0].sum() - df_alpha.iloc[:, 1].sum())
    Q_beta  = E_CHARGE * (df_beta.iloc[:, 0].sum()  - df_beta.iloc[:, 1].sum())
    T_total = cfg.n_steps * dt
    I_alpha = Q_alpha / max(T_total, 1e-15)
    I_beta  = Q_beta  / max(T_total, 1e-15)

    return {
        "dt": dt,
        "df_alpha": df_alpha,
        "df_beta":  df_beta,
        "I_alpha": I_alpha,
        "I_beta":  I_beta,
        "hist_alpha": hist_a,
        "hist_beta":  hist_b,
    }

def simulate_single(cfg, V):
    """Trace a single particle for each spin channel."""
    Gamma_pair, tau = attempt_pair_rate(
        cfg.Temperature,
        float(getattr(cfg, "activation_energy", 0.05)),
        cfg.tau_0
    )
    p_jump = float(getattr(cfg, "jump_probability", getattr(cfg, "diff_coefficient", 0.15)))
    p_jump = max(1e-4, min(0.3, p_jump))
    dt = p_jump / Gamma_pair

    df_a, hist_a, traj_a = run_channel(cfg, V, Gamma_pair, dt, spin_sigma=+1, Np=1, single_particle=True)
    df_b, hist_b, traj_b = run_channel(cfg, V, Gamma_pair, dt, spin_sigma=-1, Np=1, single_particle=True)
    return {"traj_alpha": traj_a, "traj_beta": traj_b, "dt": dt}
