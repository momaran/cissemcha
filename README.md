# CISS/eMChA Kinetic Monte Carlo (KMC) Simulator

A minimal, modular **kinetic Monte Carlo** code to simulate **spin-selective transport (CISS)** and **magneto-chiral anisotropy (eMChA)** in a 1D channel under voltage bias. The code provides a small GUI to set parameters, several ready-to-run simulations, and figures/animations as outputs.

---

## 0) Quick start

```bash
# Python 3.9–3.12
pip install numpy matplotlib seaborn imageio

# Run (opens the GUI)
python main.py

```markdown

## 1) Repository layout

cissemcha-main/
├── main.py                     # Entry point: opens the GUI and dispatches
├── read_data.py                # Tkinter GUI, validation, derived quantities
├── diffusion_mechanism.py      # Physics + KMC core (Glauber, CISS, eMChA)
├── simulations/
│   ├── Simulation_type.py      # Dispatcher (int → module.run)
│   ├── plot_style.py           # Plot style helpers (matplotlib/seaborn)
│   ├── initial_state.py        # (optional helpers for states)
│   ├── steps_counters/
│   │   └── sim.py              # Drained/Sourced vs steps (+ optional GIF over q_CISS)
│   ├── iv_three_sources/
│   │   └── sim.py              # I–V: unpolarized, alpha-only, beta-only
│   ├── iv_alpha_qciss/
│   │   └── sim.py              # I–V vs multiple q_CISS (α/β/total)
│   └── trajectory_histograms/
│       └── sim.py              # GIF: α and β position histograms along the channel
├── last_config.json            # Auto-saved GUI config (offered on next run)
└── outputs/                    # Results, time-stamped subfolders

## 2) How to run

1. **Launch**: `python main.py` (the Tk GUI appears).
2. **Choose Simulation type** (see §5).
3. **Fill parameters** (defaults provided), toggle **Enable CISS**, **Enable eMChA**, and **Self-consistent internal B(I)** if desired.
4. **Press OK**. Outputs are saved under `outputs/<timestamp>/`.

> The GUI stores your last choices in `last_config.json`. On the next run, it offers to reuse them.

---

## 3) Parameters (GUI fields)

### Core physics
- **Temperature** (K)  
- **activation_energy** \(U_\mathrm{eff}\) (eV)  
- **tau_0** attempt time (s)  
- **molecule_length** \(L\) (m)  
- **positions** number of lattice sites \(N\) (spacing \(a=L/(N-1)\))  
- **voltage_magnitude** amplitude \(V\) (V)  
- **spin_ratio** (α/β) → sets `alpha_spins` and `beta_spins`  
- **number_spins** total walkers  

### KMC numerics
- **jump_probability** \(p_\mathrm{jump}\in(0,1)\) sets \(\Delta t = p_\mathrm{jump}/\Gamma_\mathrm{pair}\).  
  (You can ignore legacy `diff_coefficient` if you use `jump_probability`.)

### CISS
- **ciss_effect** \(=q_\mathrm{CISS}\in[0,1]\) (coupling strength)  
- **chirality** \(\chi=\pm 1\)

### eMChA
- **beta_emcha** (dimensionless scale for eMChA in self-consistent mode)  
- **B_ext** external magnetic field (T)  
- **self_consistent_B** toggle; if true and `helix_pitch > 0`, add internal \(B_\text{int}=\mu_0 I/\text{pitch}\)  
- **helix_pitch** (m) for internal field  

### Advanced (per-simulation)
- **qCISS_vector** `[float, …]` for parameter sweeps in IV/steps GIF  
- **V_vector** `[float, …]` for IV scans (override default grid)  

---

## 4) Physics & model (implemented equations)

### 4.1 Lattice, field, and rates

- 1D lattice with \(N\) sites, spacing \(a=L/(N-1)\), uniform field \(E=V/L\).  
- **Per-hop field energy** (right − left):  
  \[
  \Delta E_\text{field} = eEa = e\,\frac{V}{L}\,a.
  \]  
  Implemented in `field_energy_bias_eV(...)` (in eV units).

- **Arrhenius attempt time**:  
  \[
  \tau=\tau_0 \exp\!\left(\frac{U_\text{eff}}{k_BT}\right), \quad \Gamma_\text{pair}\approx \frac{2}{\tau}.
  \]  
  Implemented in `attempt_pair_rate(...)`.

- **Glauber split** (detailed balance for hop asymmetry):  
  \[
  w_R=\frac{\Gamma}{1+e^{-\Delta E/k_BT}},\quad w_L=\Gamma-w_R.
  \]  
  Implemented in `glauber_rates(...)`.

- **KMC step**: choose \(p_\mathrm{jump}\), then \(\Delta t=p_\mathrm{jump}/\Gamma\).  
  Per step:  
  \[
  p_R=w_R\Delta t,\quad p_L=w_L\Delta t,
  \]  
  and draw right/left/no-move for each walker.  
  Probabilities are capped at <1 for stability.

#### Boundary conditions (periodic/circular)
- Right hop from \(N-1\) → counts as **Drained**, reappears at 0.  
- Left hop from 0 → counts as **Sourced**, reappears at \(N-1\).  
- **Net current** = accumulated (Drained − Sourced) / total time.  

---

### 4.2 CISS (spin-dependent splitting)

CISS is modeled as a spin-dependent energy splitting added to \(\Delta E_\text{field}\). Two parameterizations:

#### A) Exact, polarization-driven (default)
Prescribe a target spin polarization per hop \(P(V)\), e.g.  
\[
P(V)=q_\mathrm{CISS}\tanh\!\left(\eta \frac{eEa}{k_BT}\right), \quad |P|<1.
\]  
Then the energy ensuring rate asymmetry = \(P\) under Glauber is:  
\[
\Delta E_\text{CISS}(V)=2k_BT\,\operatorname{arctanh}(\chi P(V)).
\]  

- In code: `ciss_energy_split_eV(...)`  
- Near \(V=0\): **linear slope**  
- Large \(|V|\): saturation (tanh-like)  
- Coupling is `q_ciss` (GUI: *CISS Effect*).  

#### B) Linear-in-polarization (optional)
\[
\Delta E_\text{CISS}=\lambda_\text{CISS}\,k_BT\,P(V).
\]  
Equivalent to exact form at small \(P\); reproduces **gentle tilt near 0**.

#### Total bias for spin \(\sigma=\pm 1\)
\[
\Delta E_\text{total}^{(\sigma)} = \Delta E_\text{field} + \sigma\Delta E_\text{CISS} + \Delta E_\text{eMChA}.
\]

---

### 4.3 eMChA (even-in-V at fixed B)

Two models:

#### Self-consistent \(B\cdot I\) (implemented)
\[
\Delta E_\text{eMChA} = \beta_\text{emcha}\,\chi\,(B_\text{tot} I_\text{est}),
\]  
with \(B_\text{tot}=B_\text{ext}+\mu_0 I_\text{est}/\text{pitch}\) if self-consistent option is enabled.  
- Since \(I_\text{est}(V)\) flips with \(V\), \(B\cdot I\) is **even in \(V\)** — hallmark of eMChA.

#### Explicit quadratic \(BV^2\) (optional)
\[
\Delta E_\text{eMChA} = \beta^{(2)}_\text{emcha}\,\chi\,B\,V^2,
\]  
Strictly even in \(V\), convenient to separate CISS (odd) from eMChA (even).  

---

## 5) Simulations (outputs per module)

### 5.1 `steps_counters/sim.py`
- Drained/Sourced vs Steps (PNG).  
- Optional GIF sweep over `qCISS_vector`.  
- Output: `drained_vs_steps.png`, `sourced_vs_steps.png`, `drain_source_qciss.gif`.

### 5.2 `iv_three_sources/sim.py`
- I–V for **Unpolarized**, **Alpha-only**, **Beta-only**.  
- Output: `IV_three_cases.png`.

### 5.3 `iv_alpha_qciss/sim.py`
- I–V curves for multiple `qCISS_vector` (α/β/total modes).  
- Output: `IV_alpha_vs_qciss.png` etc.  

### 5.4 `trajectory_histograms/sim.py`
- GIF of α/β particle histograms over time.  
- Output: `spins_evolution.gif`.  

---

## 6) Outputs

Each run creates a folder `outputs/YYYYMMDD_HHMMSS/` with:
- Static PNGs  
- Optional GIFs  
- Optional CSV (if enabled)

---

## 7) Practical tips

- **Gentle slope near 0 V**: increase `positions`, set `q_ciss=0.05–0.2`, use `T≈300 K`, `jump_probability≈0.02–0.05`, `n_steps≥3e4`.  
- **Parity separation**:  
  - CISS (B=0): odd in V.  
  - eMChA (fixed B): even in V.  
- **Units**: energies in eV, \(k_BT\) in eV, currents in A.  
- Prefer `jump_probability`; ignore `diff_coefficient`.

---

## 8) Dispatcher

File: `simulations/Simulation_type.py`

```python
class SimulationType(Enum):
    STEPS_COUNTERS        = 0
    IV_THREE_SOURCES      = 1
    IV_ALPHA_QCISS        = 2
    TRAJECTORY_HISTOGRAMS = 3

_ROUTES = {
    0: ("simulations.steps_counters.sim",        "run"),
    1: ("simulations.iv_three_sources.sim",      "run"),
    2: ("simulations.iv_alpha_qciss.sim",        "run"),
    3: ("simulations.trajectory_histograms.sim", "run"),
}


