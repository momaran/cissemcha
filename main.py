import sys, os
from types import SimpleNamespace
from datetime import datetime
from pathlib import Path

PROJ_ROOT = Path(__file__).resolve().parent
if str(PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJ_ROOT))

from read_data import read_data
from simulations.Simulation_type import dispatch  # o 'simulation_type' si renombras

def _to_namespace(d: dict) -> SimpleNamespace:
    return SimpleNamespace(**d)

def _make_outdir(base="outputs") -> Path:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = PROJ_ROOT / base / ts
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

def main():
    cfg_dict = read_data()
    cfg = _to_namespace(cfg_dict)
    outdir = _make_outdir()
    print("[INFO] Running simulation_type:", cfg.simulation_type)
    print("[INFO] Output directory:", outdir)
    dispatch(cfg.simulation_type, cfg, outdir)
    print("[INFO] Done.")

if __name__ == "__main__":
    main()
