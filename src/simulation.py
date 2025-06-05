import numpy as np
from .config import load_config
from .diffusion import apply_diffusion_mechanism


def run_simulation(config_path="user_configurations.xlsx"):
    """Run the diffusion simulation using the provided configuration file.

    Parameters
    ----------
    config_path : str, optional
        Path to the Excel configuration file, by default "user_configurations.xlsx".

    Returns
    -------
    tuple
        summary dataframe, alpha dataframe, beta dataframe.
    """
    config = load_config(config_path)

    alpha_state = np.zeros((config["alpha_spins"], config["n_steps"]), dtype=int)
    beta_state = np.zeros((config["beta_spins"], config["n_steps"]), dtype=int)

    if config["alpha_init_position"] < config["positions"]:
        alpha_state[:, 0] = config["alpha_init_position"]
    if config["beta_init_position"] < config["positions"]:
        beta_state[:, 0] = config["beta_init_position"]

    summary_df, df_alpha, df_beta, *_ = apply_diffusion_mechanism(
        config,
        alpha_state,
        beta_state,
        config["ciss_effect"],
        config["dV"],
    )

    return summary_df, df_alpha, df_beta
