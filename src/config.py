import pandas as pd
import numpy as np

def load_config(filepath="user_configurations.xlsx"):
    """
    Loads and parses the full set of simulation parameters from an Excel file.

    The Excel file must have a header row, with the second column labeled 'Value',
    and the first column optionally describing the variable. The function assumes
    a fixed ordering of parameters as in the legacy version of the code.

    Args:
        filepath (str): Path to the user configuration Excel file.

    Returns:
        dict: Configuration dictionary with all parameters required by the simulation.
    """
    # Read only the 'Value' column, skipping the first row which may be a label
    values = pd.read_excel(filepath, usecols=[1]).iloc[1:, 0].dropna().to_numpy()

    # Extract and organize all parameters by index
    config = {
        "n_strands": int(values[0]),
        "n_turns_step": int(values[1]),
        "n_steps": int(values[2]),
        "Temperature": float(values[3]),
        "relaxation_time": float(values[4]),
        "type_of_pulse": int(values[5]),
        "voltage_magnitude": float(values[6]),
        "voltage_freq": float(values[7]),
        "n_cycles": int(values[8]),
        "alpha_init_position": int(values[9]),
        "beta_init_position": int(values[10]),
        "init_polarization_fraction": int(values[11]),
        "helix_twisting": int(values[12]),
        "magnetochiral_anisotropy": int(values[13]),
        "number_spins": int(values[14]),
        "diff_coefficient": float(values[15]),
        "Boltzmann_constant": float(values[16]),
        "electron_charge": float(values[17]),
        "ciss_effect": float(values[18]),
        "molecule_length": float(values[19]),
        "positions": int(values[20]),
        "starting_mode": float(values[21]),
        "spin_ratio": int(values[22]),
        "setting_seed": int(values[23]),
        "simulation_type": int(values[24]),
        "save_results": int(values[25])
    }

    # Derived parameters
    config["total_strands"] = config["n_strands"] * 2
    config["relaxation_rate"] = 1 / config["relaxation_time"] if config["relaxation_time"] != 0 else np.nan
    config["k"] = 11604.525  # eV/K

    # Compute number of spins for each channel
    config["alpha_spins"] = int(config["number_spins"] / (config["spin_ratio"] + 1))
    config["beta_spins"] = config["number_spins"] - config["alpha_spins"]

    # Step size and voltage gradient
    config["dx"] = config["molecule_length"] / config["positions"]
    config["dt"] = (2e-10 / (6 * config["diff_coefficient"] * config["positions"]))
    config["dV"] = config["voltage_magnitude"] / config["positions"]

    # Area estimation (arbitrary units)
    config["Area"] = config["molecule_length"] * config["positions"] * config["number_spins"]

    return config
