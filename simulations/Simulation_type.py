# simulations/Simulation_type.py
from enum import Enum
from importlib import import_module
from pathlib import Path
from typing import Tuple, Callable, Any

class SimulationType(Enum):
    STEPS_COUNTERS        = 0  # drained/sourced vs steps (+ optional GIF over q_ciss)
    IV_THREE_SOURCES      = 1  # I–V: unpolarized, beta-only, alpha-only
    IV_ALPHA_QCISS        = 2  # I–V: alpha-only for several q_ciss
    TRAJECTORY_HISTOGRAMS = 3  # GIF: particle histograms along the channel (α & β)
    IV_EMCHA             = 4  # I–V: with electrical magnetochiral anisotropy (EMCA)

# Map enum -> (module_path, function_name)
_ROUTES = {
    SimulationType.STEPS_COUNTERS.value:        ("simulations.steps_counters.sim",        "run"),
    SimulationType.IV_THREE_SOURCES.value:      ("simulations.iv_three_sources.sim",      "run"),
    SimulationType.IV_ALPHA_QCISS.value:        ("simulations.iv_alpha_qciss.sim",        "run"),
    SimulationType.TRAJECTORY_HISTOGRAMS.value: ("simulations.trajectory_histograms.sim", "run"),
    SimulationType.IV_EMCHA.value:               ("simulations.iv_emcha.sim",               "run"),
}

def _resolve_route(simulation_type: int) -> Tuple[str, str]:
    if simulation_type not in _ROUTES:
        valid = ", ".join(str(k) for k in sorted(_ROUTES.keys()))
        raise ValueError(f"Unknown simulation_type={simulation_type}. Valid keys: {valid}")
    return _ROUTES[simulation_type]

def dispatch(simulation_type: int, cfg: Any, outdir) -> Any:
    """
    Dynamically import the selected simulation and run it.
    Ensures 'outdir' is a pathlib.Path so downstream code can safely call mkdir().
    """
    module_path, fn_name = _resolve_route(simulation_type)

    # Import the module and fetch the run() function
    mod = import_module(module_path)
    fn: Callable = getattr(mod, fn_name)

    # Normalize outdir to Path (robust to str being passed in)
    outdir = Path(outdir)

    # Create the directory upfront (ok if it already exists)
    outdir.mkdir(parents=True, exist_ok=True)

    # Call the simulation entry point
    return fn(cfg, outdir)
