import pandas as pd

# Map numeric simulation types to readable mode names
SIMULATION_MODES = {
    1: "polarization_vs_voltage",
    2: "polarization_vs_position",
    3: "spin_density_evolution",
    # Add more modes as needed
}

def load_config(filepath="user_configurations.xlsx"):
    """
    Load simulation configuration parameters from an Excel file.

    The Excel file must have two columns: 'Variable' and 'Value'.

    Args:
        filepath (str): Path to the Excel file with simulation parameters.

    Returns:
        dict: A dictionary containing simulation parameters.
              Includes a key 'mode' with a readable name.
    """
    df = pd.read_excel(filepath)
    config = {}

    # Parse each row into the config dictionary
    for _, row in df.iterrows():
        key = row["Variable"]
        value = row["Value"]
        config[key] = value

    # Convert simulation_type (int) to a readable mode string
    sim_type = int(config.get("simulation_type", 1))
    config["mode"] = SIMULATION_MODES.get(sim_type, "unknown_mode")

    return config
