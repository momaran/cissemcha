import numpy as np
import pandas as pd

def diffusion_mechanism(spin_type, config, state_matrix, ciss_effect, dV):
    """
    Simulates spin-dependent diffusion for either α (spin -1) or β (spin +1) particles.

    Args:
        spin_type (int): -1 for α spins, +1 for β spins.
        config (Namespace): Configuration object containing simulation parameters.
        state_matrix (ndarray): (n_spins, n_steps) matrix of particle positions.
        ciss_effect (float): Magnitude of the CISS effect.
        dV (float): Voltage difference per step (longitudinal voltage gradient).

    Returns:
        tuple:
            - updated state_matrix (ndarray)
            - total drained spins (int)
            - total sourced spins (int)
            - dataframe of spin tracking per step (pd.DataFrame)
            - average right probability (float)
            - average left probability (float)
            - net spin current (float)
    """
    if config.setting_seed == 1:
        np.random.seed(5)

    n_steps = config.n_steps
    positions = config.positions
    diff_coeff = config.diff_coefficient
    max_position = positions

    drained_spins = np.zeros(n_steps, dtype=int)
    sourced_spins = np.zeros(n_steps, dtype=int)

    # OPTION: tanh-based spin-dependent modulation
    ciss_contribution = ciss_effect * spin_type * np.tanh(10 * dV)

    right_prob = (1 / (1 + np.exp(-2 * 5 * dV * (1 + ciss_contribution)))) * diff_coeff
    left_prob = diff_coeff - right_prob

    for t in range(1, n_steps):
        random_values = np.random.rand(len(state_matrix))
        prev = state_matrix[:, t - 1]

        # Right moves
        right_mask = random_values <= right_prob
        move_right = right_mask & (prev != max_position)
        wrap_right = right_mask & (prev == max_position)

        # Left moves
        left_mask = (random_values > right_prob) & (random_values <= diff_coeff)
        move_left = left_mask & (prev != 0)
        wrap_left = left_mask & (prev == 0)

        # Update positions
        state_matrix[move_right, t] = prev[move_right] + 1
        state_matrix[wrap_right, t] = 0
        drained_spins[t] += np.sum(wrap_right)

        state_matrix[move_left, t] = prev[move_left] - 1
        state_matrix[wrap_left, t] = max_position
        sourced_spins[t] += np.sum(wrap_left)

        # Unmoved particles
        stay = ~(move_right | wrap_right | move_left | wrap_left)
        state_matrix[stay, t] = prev[stay]

    # Compute spin currents (in arbitrary units)
    I_drain = drained_spins.sum() / n_steps
    I_source = sourced_spins.sum() / n_steps
    I_difference = I_drain - I_source

    spin_label = "α" if spin_type == -1 else "β"
    df = pd.DataFrame({
        "Step": np.arange(n_steps),
        f"Drained {spin_label} spins": drained_spins,
        f"Sourced {spin_label} spins": sourced_spins,
        f"Right Probability ({spin_label})": right_prob,
        f"Left Probability ({spin_label})": left_prob
    })

    return state_matrix, drained_spins.sum(), sourced_spins.sum(), df, right_prob, left_prob, I_difference


def apply_diffusion_mechanism(config, alpha_state_matrix, beta_state_matrix, ciss_effect, dV):
    """
    Runs the spin diffusion simulation for both α and β spins and aggregates results.

    Args:
        config (Namespace): Configuration with all parameters.
        alpha_state_matrix (ndarray): Initial positions for α spins.
        beta_state_matrix (ndarray): Initial positions for β spins.
        ciss_effect (float): Value of CISS parameter.
        dV (float): Voltage step per site.

    Returns:
        tuple:
            - summary dataframe (pd.DataFrame)
            - dataframe for α spins
            - dataframe for β spins
            - updated alpha_state_matrix
            - updated beta_state_matrix
    """
    # α spins (spin_type = -1)
    alpha_state_matrix, drained_α, sourced_α, df_α, r_prob_α, l_prob_α, I_α = diffusion_mechanism(
        -1, config, alpha_state_matrix, ciss_effect, dV
    )

    # β spins (spin_type = +1)
    beta_state_matrix, drained_β, sourced_β, df_β, r_prob_β, l_prob_β, I_β = diffusion_mechanism(
        +1, config, beta_state_matrix, ciss_effect, dV
    )

    summary_df = pd.DataFrame({
        "Spin Type": ["α", "β"],
        "Total Drained": [drained_α, drained_β],
        "Total Sourced": [sourced_α, sourced_β],
        "Avg Right Probability": [np.mean(r_prob_α), np.mean(r_prob_β)],
        "Avg Left Probability": [np.mean(l_prob_α), np.mean(l_prob_β)],
        "Net Current": [I_α, I_β]
    })

    return summary_df, df_α, df_β, alpha_state_matrix, beta_state_matrix
