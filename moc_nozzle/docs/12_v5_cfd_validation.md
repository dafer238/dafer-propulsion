# V5 — CFD Validation

## Why Validate Against CFD?

The MOC solver makes several assumptions that simplify the physics significantly:

1. **Inviscid flow** — no viscosity, no friction (partially corrected by the V4 displacement-body
   model, but that correction is itself an approximation)
2. **Isentropic flow** — no shocks, no entropy generation; the method of characteristics breaks down
   at any shock
3. **Frozen or equilibrium chemistry** — species concentrations are fixed or assumed to equilibrate
   instantaneously (V3 addresses the thermodynamic properties but not the kinetics)
4. **Calorically perfect gas** — constant γ everywhere in the flow (V2 partially fixes this, but
   the MOC kernel still uses a single γ)
5. **Steady, axisymmetric flow** — no injector asymmetry, no transient startup effects

CFD (Computational Fluid Dynamics) with a Navier-Stokes solver relaxes assumptions 1 and 2
simultaneously and can run in a fully 2D axisymmetric or 3D domain. Solving the Reynolds-Averaged
Navier-Stokes (RANS) equations with a turbulence model adds viscous effects beyond the simple
displacement-body approximation of V4.

Running a CFD simulation on your MOC-designed nozzle contour will reveal:

- Whether the MOC contour **actually produces shock-free flow** at the design Mach — an incorrect
  contour will show oblique shocks emanating from the throat region
- How large the **viscous losses** are and how well the V4 estimate compares to the full viscous
  simulation
- Whether the **design exit Mach number is achieved** across the full exit plane (not just on the
  axis)
- **Near-wall flow behavior**: wall heating rates, risk of flow separation in the diverging section,
  and boundary layer growth rates
- The accuracy of the **thrust coefficient C_F** and **specific impulse Isp** predictions

CFD validation is the gold standard for nozzle design verification before hardware is built. A
computation that costs a few hours on a laptop can catch a design error that would cost weeks and
thousands of dollars to discover during hot-fire testing.

---

## Axisymmetric 2D CFD vs. Full 3D

For a circularly symmetric nozzle (the type produced by this MOC code), you can run a
**2D axisymmetric** CFD simulation instead of a full 3D one:

- The computational domain is a 2D mesh in the r-x plane with an axisymmetric boundary condition
  on the centerline (r = 0)
- The solver internally applies the geometric source terms that account for the cylindrical
  coordinate system
- This is **10–100× cheaper** than a full 3D simulation in terms of both mesh size and
  computational time
- For a simple axisymmetric nozzle, the result is **identical** to 3D (no information is lost)

The axisymmetric domain looks like this in cross-section:

```
  r
  ↑
  |  [chamber]──[convergent]──[throat]──[divergent MOC contour]──[exit]  ← wall BC
  |
  |   (flow domain: subsonic chamber + sonic throat + supersonic divergent)
  |
  0─────────────────────────────────────────────────────────────────────► x
     axis of symmetry (r = 0): symmetry / axisymmetric BC
```

**When to use full 3D:**
- Injector-head asymmetry (off-axis fuel/oxidizer injection patterns causing non-uniform chamber flow)
- Thrust vector control (TVC) deflections creating geometric asymmetry
- Non-circular cross-sections (plug nozzles, rectangular aerospike concepts)
- Combustion instability coupling studies (require 3D + acoustic modes)

For all standard nozzle design validation in this project, 2D axisymmetric is appropriate.

---

## Recommended Free CFD Tools

### SU2 (Stanford University Unstructured)

- **Website:** https://su2code.github.io
- **Language:** C++ solver, Python scripting interface
- **License:** LGPL-2.1 (free, open-source)
- **Strengths:** Euler (inviscid) and RANS (viscous turbulent) solvers, adjoint-based shape
  optimization, excellent documentation, built-in axisymmetric mode, active development community
- **Input format:** `.su2` unstructured mesh file + `.cfg` plain-text configuration file
- **Output:** `.vtu` / `.csv` flow field files readable in ParaView
- **Recommended for this project:** **best overall fit** for nozzle validation — it directly
  supports inviscid Euler for contour verification, then viscous RANS to quantify the boundary layer
- **Getting started:** https://su2code.github.io/docs/home/

### OpenFOAM

- **Website:** https://www.openfoam.com (ESI release) or https://openfoam.org (Foundation release)
- **Language:** C++ solver, case-directory structure with dictionary files
- **License:** GPL-3.0 (free, open-source)
- **Strengths:** enormous user community, large library of turbulence models and thermophysical
  property databases, mature mesh motion and multi-phase capabilities
- **Relevant solvers:** `sonicFoam` (transonic/supersonic density-based), `rhoCentralFoam`
  (central-scheme compressible), `rhoSimpleFoam` (RANS steady state)
- **Axisymmetric:** use a "wedge" mesh (a thin angular slice of the domain, typically 5°) with
  `wedge` boundary conditions on the two azimuthal faces; the solver treats it as axisymmetric
- **Learning curve:** steeper than SU2 due to the dictionary-based case setup, but extremely
  flexible once learned

### ANSYS Fluent / CFX

- **License:** Commercial (expensive for professional use); a free Student license is available
  at https://www.ansys.com/academic/students with limited mesh size (~512k cells)
- **Strengths:** industry-standard, excellent GUI, very mature density-based compressible solver in
  Fluent with proven turbulence models
- **Axisymmetric:** built into Fluent; select "axisymmetric" in the 2D space option with a
  density-based solver and ideal-gas equation of state
- **Recommended:** if you already have access through an institution, Fluent is easy to set up for
  nozzle flows; otherwise SU2 or OpenFOAM are preferred for open-source work

### CONVERGE CFD

- Commercial, high cost, primarily used for internal combustion engine simulations. Not recommended
  for this project.

---

## Setting Up the CFD Simulation

### Step 1: Export the Nozzle Contour from Rust

Add an export function to `src/geometry/wall.rs` or directly in `main.rs`:

```rust
/// Export wall contour points to a CSV file for use by GMSH and other tools.
/// Points should be normalized to r_t = 1, or in physical meters — document the convention.
///
/// # Arguments
/// * `path`   - Output file path (e.g., "output/wall_divergent.csv")
/// * `points` - Slice of (x, r) pairs from the MOC wall nodes
pub fn export_wall_csv(path: &str, points: &[(f64, f64)]) {
    use std::io::Write;
    let mut f = std::fs::File::create(path)
        .expect("Cannot create wall CSV file");
    writeln!(f, "x,r").unwrap();
    for (x, r) in points {
        writeln!(f, "{:.8},{:.8}", x, r).unwrap();
    }
    println!("Exported {} wall points to {}", points.len(), path);
}
```

Call this in `main.rs` after computing the wall:

```rust
let wall_pts: Vec<(f64, f64)> = nozzle.solver
    .wall_nodes()
    .iter()
    .map(|n| (n.x, n.y))
    .collect();

export_wall_csv("output/wall_divergent.csv", &wall_pts);
```

The resulting `wall_divergent.csv` has two columns `x` and `r` (both in meters if you scaled your
MOC by r_t). This file is the primary input to the mesh generator.

### Step 2: Define the Full Nozzle Geometry

The CFD mesh requires the complete nozzle profile from the chamber to the exit, not just the
divergent section. The geometry consists of the following curve segments (listed clockwise from the
inlet):

1. **Inlet face** — a vertical line from the axis to the wall at x = x_inlet (the chamber entrance
   plane)
2. **Chamber wall** — a short horizontal or slightly converging wall at r = r_chamber
3. **Convergent section** — from r_chamber to r_throat (from V6 `convergent.rs`, or use a circular
   arc + straight section)
4. **Throat radius region** — a circular arc of radius R_c (typically 0.5–1.5 × r_t) blending the
   convergent into the divergent; this ensures a smooth sonic line
5. **Divergent MOC wall** — the spline through the exported `wall_divergent.csv` points
6. **Exit face** — a vertical line from the axis to the wall at x = x_exit
7. **Axis** — the centerline from inlet to exit at r = 0

```
r
↑     [2] chamber      [3] convergent   [4] throat   [5] divergent MOC
│  ─────────────────────────────────────────────────────────────────────  ← wall
│                                           ╭─────╮
│                                     ╭────╯       ╰────────────────────
│ [1]                                                               [6]
│ inlet                                                             exit
│
0──────────────────────────────────────────────────────────────────────► x
  [7] axis (r = 0)
```

### Step 3: Mesh Generation with GMSH

GMSH (https://gmsh.info) is the recommended free meshing tool. It can import curve data from CSV,
generate structured or unstructured triangular / quadrilateral meshes, and export directly to the
`.su2` format required by SU2.

Key meshing guidelines for a nozzle flow:

| Region | Recommended cell size | Rationale |
|--------|----------------------|-----------|
| Throat (±1 r_t) | r_t / 100 to r_t / 50 | Highest gradients; sonic line |
| Wall boundary layer | first cell height ≈ y+ ~1 for RANS | Resolves viscous sublayer |
| Near-wall growth rate | 1.1–1.15 per layer | Smooth expansion into freestream |
| Freestream core | r_t / 10 | Lower gradients; save cells |
| Exit plane | r_t / 20 | Need good resolution to extract exit profile |

A mesh of ~10,000–50,000 cells is typically sufficient for a 2D axisymmetric nozzle validation.

See the GMSH Python script template section below for a starting point.

### Step 4: Configure SU2 for Euler (Inviscid) Flow

Start with an **Euler** (inviscid) simulation to validate the MOC contour shape itself, before
adding the complexity of viscous effects. If the Euler CFD shows shocks, the MOC contour is wrong;
add viscosity only once the inviscid solution is clean.

Save the following as `nozzle_euler.cfg`:

```
% ============================================================
% SU2 Configuration — Euler (inviscid) nozzle validation
% ============================================================

% --- Physics ---
SOLVER                  = EULER
KIND_TURB_MODEL         = NONE
RESTART_SOL             = NO

% --- Thermodynamic model ---
FLUID_MODEL             = IDEAL_GAS
GAMMA_VALUE             = 1.20       % Use γ from CEA / V3 output
GAS_CONSTANT            = 692.0      % R_specific [J/(kg·K)] for LOX/LH₂ mixture

% --- Axisymmetric ---
AXISYMMETRIC            = YES

% --- Boundary conditions ---
% Inlet: specify stagnation (total) temperature and pressure
MARKER_INLET            = ( inlet, 3500.0, 10000000.0 )   % T0 [K], P0 [Pa]

% Exit: supersonic outflow — do not specify anything, the flow sets its own state
% (or specify a back pressure for subsonic exit, which is not the case here)
MARKER_OUTLET           = ( outlet, 60000.0 )   % static pressure [Pa]; only active
                                                 % if flow becomes subsonic at exit

% Symmetry axis
MARKER_SYM              = ( axis )

% Nozzle wall: slip (inviscid) wall
MARKER_EULER            = ( wall )

% --- Numerical scheme ---
CONV_NUM_METHOD_FLOW    = ROE         % Roe upwind scheme — good for shocks
MUSCL_FLOW              = YES         % 2nd-order reconstruction
SLOPE_LIMITER_FLOW      = VENKATAKRISHNAN  % prevents spurious oscillations near shocks
CFL_NUMBER              = 1.0
LINEAR_SOLVER           = FGMRES
LINEAR_SOLVER_ERROR     = 1E-10
MAX_ITER                = 5000

% --- Output ---
OUTPUT_FILES            = ( RESTART, PARAVIEW, SURFACE_CSV )
CONV_FIELD              = RMS_DENSITY
CONV_RESIDUAL_MINVAL    = -10        % converge density residual to 10^{-10}
OUTPUT_WRT_FREQ         = 250
VOLUME_OUTPUT           = ( MACH, PRESSURE, TEMPERATURE, DENSITY, VELOCITY )
SURFACE_OUTPUT          = ( MACH, PRESSURE, TEMPERATURE )
```

**Running SU2 Euler:**

```sh
# Single core
SU2_CFD nozzle_euler.cfg

# Parallel (4 cores)
mpirun -n 4 SU2_CFD nozzle_euler.cfg
```

The output `flow.vtu` can be opened in ParaView. Look for:
- Mach number contour: should smoothly increase from ~0 in the chamber to M_exit at the exit
- Pressure contour: should smoothly decrease with no discontinuities (shocks appear as sharp lines)
- If you see oblique shocks in the divergent section, the MOC contour needs correction

### Step 5: Configure SU2 for RANS (Viscous) Flow

Once the Euler simulation confirms the contour is shock-free, run a viscous RANS simulation to
quantify boundary layer effects and compare with the V4 estimate.

Save the following as `nozzle_rans.cfg` (differences from the Euler config are marked):

```
% ============================================================
% SU2 Configuration — RANS (viscous) nozzle validation
% ============================================================

% --- Physics (CHANGED from EULER) ---
SOLVER                  = RANS
KIND_TURB_MODEL         = SST        % Menter SST k-ω — best for nozzle flows
                                     % (good near-wall and freestream behavior)
RESTART_SOL             = NO

% --- Thermodynamic model ---
FLUID_MODEL             = IDEAL_GAS
GAMMA_VALUE             = 1.20
GAS_CONSTANT            = 692.0

% --- Viscosity model (ADDED for RANS) ---
VISCOSITY_MODEL         = SUTHERLAND
MU_REF                  = 1.0E-5    % reference viscosity [Pa·s]
MU_T_REF                = 1000.0    % reference temperature [K]
SUTHERLAND_CONSTANT     = 200.0     % Sutherland constant S [K]

% --- Axisymmetric ---
AXISYMMETRIC            = YES

% --- Boundary conditions ---
MARKER_INLET            = ( inlet, 3500.0, 10000000.0 )
MARKER_OUTLET           = ( outlet, 60000.0 )
MARKER_SYM              = ( axis )

% Wall: no-slip (viscous) wall, isothermal at T_wall = 700 K
% (for adiabatic wall use MARKER_HEATFLUX instead with heat flux = 0)
MARKER_ISOTHERMAL       = ( wall, 700.0 )   % T_wall [K]

% --- Turbulence inlet conditions ---
% Turbulent intensity 5%, viscosity ratio μ_t/μ = 10 (typical for rocket chambers)
FREESTREAM_TURBULENCEINTENSITY = 0.05
FREESTREAM_TURB2LAMVISCRATIOITY = 10.0

% --- Numerical scheme (same as Euler) ---
CONV_NUM_METHOD_FLOW    = ROE
MUSCL_FLOW              = YES
SLOPE_LIMITER_FLOW      = VENKATAKRISHNAN
CFL_NUMBER              = 0.5        % lower than Euler; viscous terms need smaller CFL
LINEAR_SOLVER           = FGMRES
MAX_ITER                = 10000      % viscous simulations need more iterations

% --- Output ---
OUTPUT_FILES            = ( RESTART, PARAVIEW, SURFACE_CSV )
CONV_FIELD              = RMS_DENSITY
CONV_RESIDUAL_MINVAL    = -8
OUTPUT_WRT_FREQ         = 500
VOLUME_OUTPUT           = ( MACH, PRESSURE, TEMPERATURE, DENSITY,
                             VELOCITY, LAMINAR_VISCOSITY, EDDY_VISCOSITY )
SURFACE_OUTPUT          = ( MACH, PRESSURE, TEMPERATURE, SKIN_FRICTION )
```

**Mesh requirement for RANS:** The first cell height off the wall must be chosen so that y+ ≈ 1
for the SST model (which resolves the viscous sublayer). A rough estimate for first cell height:

```
y_1 ≈ y+ · μ / (ρ · u_τ)
u_τ = √(τ_w / ρ)  where τ_w ≈ C_f · (1/2) · ρ · u²
```

For the nozzle wall at M = 2, P = 1 MPa, T = 1500 K:
- u_e ≈ 1800 m/s, ρ ≈ 2.3 kg/m³, C_f ≈ 0.002
- τ_w ≈ 0.002 × 0.5 × 2.3 × 1800² ≈ 7500 Pa
- u_τ ≈ √(7500/2.3) ≈ 57 m/s
- μ ≈ 6 × 10⁻⁵ Pa·s (at 1500 K)
- y_1 = 1 × 6×10⁻⁵ / (2.3 × 57) ≈ 4.6 × 10⁻⁷ m ≈ 0.46 μm

This is very small. Use at least 20 boundary layer cells with a growth ratio of 1.1–1.2 to
capture the full boundary layer profile.

---

## Key Metrics to Compare

After running CFD, extract and compare these quantities against the MOC predictions:

| Metric | MOC prediction | How to extract from CFD | Acceptable delta |
|--------|---------------|------------------------|-----------------|
| Exit Mach number M_exit | From area ratio / isentropic table | Area-averaged M at exit plane in ParaView | < 1% |
| Exit flow uniformity | Uniform (M = const across exit) | Standard deviation of M at exit | < 2% of M_exit |
| Shock presence | None (inviscid design) | Schlieren-style ∣∇ρ∣ plot; density gradient contours | No shock |
| Thrust coefficient C_F | Isentropic C_F formula | Integrate wall pressure + momentum flux at exit | < 2% |
| Wall pressure P_wall(x) | Isentropic P(M(x)) | Extract wall boundary from SU2 surface output | ~5% in BL region |
| Displacement thickness δ* | V4 estimate | Integrate velocity profile normal to wall | Within 30% |
| Effective exit angle | 0° (parallel flow by MOC design) | Area-weighted mean flow angle at exit | < 1° |

### Extracting Exit Mach Profile in ParaView

1. Open `flow.vtu` in ParaView
2. Add a **Plot Over Line** filter: set the line from (x_exit, 0) to (x_exit, r_exit) — along the
   exit plane
3. In the filter settings, select "Mach" as the field variable
4. Export the profile as CSV and compare to the design M_exit

**Computing area-averaged exit Mach:**

In ParaView, use **Integrate Variables** on the exit boundary extracted with **Extract Block** or
**Slice**: the result gives total area and integrated Mach, from which the mean follows.

Alternatively, read the SU2 surface CSV output: `surface_flow.csv` contains one row per boundary
node on each marked surface. Filter by `inlet` or `outlet` and compute the mean.

### Extracting Wall Pressure

SU2 writes `surface_flow.csv` with pressure at each wall boundary node. A Python snippet to
compare with the isentropic MOC prediction:

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load CFD wall pressure
cfd = pd.read_csv("surface_flow.csv")
wall = cfd[cfd["Marker"] == "wall"].sort_values("x")

# Load MOC isentropic prediction
moc = pd.read_csv("output/wall_isentropic.csv")  # x, M, P exported from Rust

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(moc["x"], moc["P"] / 1e6,  label="MOC isentropic")
axes[0].plot(wall["x"], wall["Pressure"] / 1e6, label="CFD RANS")
axes[0].set_xlabel("x [m]")
axes[0].set_ylabel("Wall static pressure [MPa]")
axes[0].legend()

axes[1].plot(moc["x"], moc["M"],         label="MOC")
axes[1].plot(wall["x"], wall["Mach"],    label="CFD RANS")
axes[1].set_xlabel("x [m]")
axes[1].set_ylabel("Mach number")
axes[1].legend()

plt.tight_layout()
plt.savefig("cfd_vs_moc_comparison.pdf")
```

---

## What to Look For: Failure Modes

Understanding what can go wrong in the simulation — and how each failure mode appears in CFD output
— is essential for diagnosing whether the issue is in the MOC design or in the CFD setup.

### Shock Waves Inside the Nozzle

**What it looks like in CFD:** Sharp discontinuities (lines) in the Mach number contour plot.
The density gradient magnitude (`||∇ρ||`) shows bright streaks — this is the "numerical schlieren"
image and directly mimics what optical schlieren diagnostics show in experimental nozzle testing.

**Cause:** The nozzle wall contour does not correctly cancel the left-running characteristics at the
wall. This can happen if:
- The design Mach number M_exit is inconsistent with the area ratio used in MOC
- The initial data line (the throat region starting conditions) uses an incorrect γ or M distribution
- A coding bug in the wall angle calculation or characteristic intersection

**Fix:**
1. Re-examine the initial data line and confirm the Mach distribution is correct for the given γ
2. Check that the wall angle at the throat (θ_max) is within the valid range: θ_max ≤ Prandtl-Meyer
   angle difference between M=1 and M_exit
3. Increase n_chars (number of characteristic lines) to get a finer approximation of the
   continuous wall

### Over-Expansion at the Exit

**What it looks like:** The exit Mach number from CFD is higher than the design M_exit, and
the flow re-compresses outside the nozzle through oblique shocks visible in the plume.

**Cause:** The effective area ratio is larger than expected. This can occur if:
- The wrong γ was used to compute A_e/A_t → M_exit (use the CEA γ from V3, not γ = 1.4)
- The boundary layer displacement (V4) was not accounted for in the hardware wall design
- The convergent section contributes extra area not accounted for

**Fix:** Recompute A_e/A_t from the correct γ. If using V4 corrections, verify that the physical
wall is `y_hardware = y_MOC + δ*(x)`, not just `y_MOC`.

### Boundary Layer Separation

**What it looks like in CFD:** In a RANS simulation, the near-wall axial velocity becomes
negative (flow reversal / recirculation) in some region of the divergent section. Streamlines near
the wall turn around. The skin friction coefficient C_f goes negative.

**Cause:** The adverse pressure gradient in the diverging section is too severe for the boundary
layer to remain attached. This typically occurs when:
- The wall divergence half-angle exceeds ~20–25° anywhere along the divergent section
- The θ_max (maximum wall angle at the throat) is too large
- The nozzle length is too short for the given area ratio (area ratio achieved in fewer
  characteristic lines means steeper wall angles)

**Fix:**
- Reduce θ_max in the MOC input (use ~15° for a mild contour, ~20° for a compact design)
- Increase n_chars to get a longer, more gradual wall contour
- Check that the wall curvature is smooth everywhere (no kinks in the spline fit)

### Exit Flow Non-Uniformity

**What it looks like:** In the exit-plane Mach profile, M varies significantly from axis to wall
(e.g., M = 3.0 at the axis vs. M = 2.8 near the wall).

**Cause:**
- Too few characteristic lines: the MOC wall is a coarse polygonal approximation to the smooth
  ideal wall, and the discontinuities in wall angle generate small waves that reach the exit plane
- Numerical diffusion in the CFD mesh if the mesh is too coarse near the exit
- Real viscous effects: even a perfectly designed inviscid contour will have slight non-uniformity
  in a viscous simulation due to the boundary layer profile

**Fix:**
- Increase n_chars (try 50 → 100 → 200 and monitor convergence of exit uniformity)
- Refine the CFD mesh near the exit plane
- Accept 1–2% non-uniformity as inherent to viscous flow and consistent with V4 corrections

### Numerical Divergence or Non-Convergence

**What it looks like:** SU2 residuals do not decrease, or the simulation crashes with NaN values.

**Cause and fix:**
- CFL too high: reduce CFL from 1.0 to 0.5 or 0.25 and restart
- Initial conditions far from the solution: use a Mach-number ramp (start the simulation
  with a lower inlet Mach and gradually increase) or run the Euler solution first and restart from it
- Mesh quality issues: check the mesh in GMSH for negative Jacobians (inverted cells), especially
  near the throat
- Supersonic outlet BC: if the outlet is supersonic, use `MARKER_OUTLET` with a very low back
  pressure (< ambient); do not over-constrain a supersonic exit

---

## Validation Workflow as a Rust Post-Processing Tool

It is useful to have the Rust binary export everything a CFD setup needs in one go: the divergent
wall, the convergent wall, and a combined full-profile file. Add this function to
`src/geometry/wall.rs` or a new `src/output/cfd_export.rs` module:

```rust
/// Export nozzle geometry for CFD meshing.
///
/// Writes three CSV files:
///   - `{prefix}_wall_convergent.csv`  — upstream convergent section
///   - `{prefix}_wall_divergent.csv`   — downstream MOC divergent contour
///   - `{prefix}_wall_full.csv`        — complete profile, chamber to exit
///
/// All coordinates in meters. x is the axial direction; r is the radial direction.
///
/// # Arguments
/// * `conv_pts` - (x, r) pairs for the convergent section, throat-to-chamber direction
/// * `div_pts`  - (x, r) pairs for the divergent MOC wall, throat to exit
/// * `prefix`   - Output file prefix, e.g. "output/nozzle"
pub fn export_for_cfd(
    conv_pts: &[(f64, f64)],
    div_pts:  &[(f64, f64)],
    prefix:   &str,
) {
    let write_csv = |path: &str, pts: &[(f64, f64)]| {
        use std::io::Write;
        let mut f = std::fs::File::create(path)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", path, e));
        writeln!(f, "x,r").unwrap();
        for (x, r) in pts {
            writeln!(f, "{:.8},{:.8}", x, r).unwrap();
        }
        println!("  Wrote {} points to {}", pts.len(), path);
    };

    println!("Exporting CFD geometry files with prefix '{}':", prefix);

    write_csv(&format!("{}_wall_convergent.csv", prefix), conv_pts);
    write_csv(&format!("{}_wall_divergent.csv",  prefix), div_pts);

    // Combine: sort by x to ensure a continuous wall from inlet to exit.
    // The convergent section typically has negative or zero x (throat at x=0),
    // and the divergent section has positive x.
    let mut combined: Vec<(f64, f64)> = conv_pts.to_vec();
    combined.extend_from_slice(div_pts);
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    combined.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12);  // remove duplicate throat point

    write_csv(&format!("{}_wall_full.csv", prefix), &combined);
}
```

Call this in `main.rs` after computing both sections:

```rust
geometry::wall::export_for_cfd(
    &convergent_wall_pts,
    &divergent_wall_pts,
    "output/nozzle",
);
// Produces: output/nozzle_wall_convergent.csv
//           output/nozzle_wall_divergent.csv
//           output/nozzle_wall_full.csv
```

---

## GMSH Python Script Template

The following Python script template reads `nozzle_wall_full.csv` and produces a 2D axisymmetric
mesh ready for SU2. You will need the `gmsh` Python package (`pip install gmsh`) and a GMSH
installation (https://gmsh.info/#Download).

This is provided as a template to adapt — the GMSH Python API is well-documented at
https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-API.

```python
"""
gmsh_nozzle.py — Generate a 2D axisymmetric nozzle mesh for SU2.

Usage:
    python gmsh_nozzle.py

Outputs:
    nozzle.su2   — unstructured mesh for SU2
    nozzle.msh   — GMSH native format (for inspection in GMSH GUI)

Requirements:
    pip install gmsh pandas numpy
"""

import gmsh
import pandas as pd
import numpy as np

# ── User settings ───────────────────────────────────────────────────────────
WALL_CSV     = "output/nozzle_wall_full.csv"  # from Rust export_for_cfd()
MESH_SIZE_THROAT  = 0.0005    # [m] cell size at the throat — tightest region
MESH_SIZE_WALL    = 0.001     # [m] cell size along the wall (away from throat)
MESH_SIZE_AXIS    = 0.002     # [m] cell size on the axis
MESH_SIZE_EXIT    = 0.001     # [m] cell size at the exit plane
BL_LAYERS         = 20        # number of boundary layer (structured quad) layers
BL_RATIO          = 1.15      # growth ratio between BL layers
BL_FIRST_HEIGHT   = 5e-7      # [m] first cell height off wall (for y+ ≈ 1 at M~2)
# ────────────────────────────────────────────────────────────────────────────

def main():
    gmsh.initialize()
    gmsh.model.add("nozzle")

    # ── Load wall profile ──────────────────────────────────────────────────
    df = pd.read_csv(WALL_CSV)
    xs = df["x"].values
    rs = df["r"].values

    # Add wall points to GMSH geometry
    wall_pts = []
    throat_idx = np.argmin(rs)   # throat is at minimum r

    for i, (x, r) in enumerate(zip(xs, rs)):
        # Mesh size: finer at the throat region
        dist_from_throat = abs(x - xs[throat_idx])
        lc = MESH_SIZE_THROAT + (MESH_SIZE_WALL - MESH_SIZE_THROAT) * min(
            dist_from_throat / (0.005), 1.0
        )
        pt = gmsh.model.geo.addPoint(x, r, 0.0, lc)
        wall_pts.append(pt)

    # Axis points: inlet and exit on the centerline
    x_inlet = xs[0]
    x_exit  = xs[-1]
    r_inlet = rs[0]    # inlet wall radius (chamber)

    pt_axis_inlet = gmsh.model.geo.addPoint(x_inlet, 0.0, 0.0, MESH_SIZE_AXIS)
    pt_axis_exit  = gmsh.model.geo.addPoint(x_exit,  0.0, 0.0, MESH_SIZE_AXIS)

    # ── Curves ────────────────────────────────────────────────────────────
    # Wall spline (all wall points, inlet to exit)
    wall_spline = gmsh.model.geo.addSpline(wall_pts)

    # Exit face: from wall exit point to axis exit point
    exit_line   = gmsh.model.geo.addLine(wall_pts[-1], pt_axis_exit)

    # Axis: from axis exit point back to axis inlet point
    axis_line   = gmsh.model.geo.addLine(pt_axis_exit, pt_axis_inlet)

    # Inlet face: from axis inlet point up to wall inlet point
    inlet_line  = gmsh.model.geo.addLine(pt_axis_inlet, wall_pts[0])

    # ── Surface ───────────────────────────────────────────────────────────
    curve_loop = gmsh.model.geo.addCurveLoop(
        [wall_spline, exit_line, axis_line, inlet_line]
    )
    surface = gmsh.model.geo.addPlaneSurface([curve_loop])

    gmsh.model.geo.synchronize()

    # ── Boundary layer (quads off wall) ──────────────────────────────────
    # GMSH BoundaryLayer field creates structured quad layers off the wall.
    f = gmsh.model.mesh.field
    bl_field = f.add("BoundaryLayer")
    f.setNumbers(bl_field, "CurvesList", [wall_spline])
    f.setNumber(bl_field,  "Size",        BL_FIRST_HEIGHT)
    f.setNumber(bl_field,  "Ratio",       BL_RATIO)
    f.setNumber(bl_field,  "NbLayers",    BL_LAYERS)
    f.setNumber(bl_field,  "Quads",       1)
    gmsh.model.mesh.field.setAsBackgroundMesh(bl_field)

    # ── Physical groups (boundary labels for SU2) ─────────────────────────
    gmsh.model.addPhysicalGroup(1, [wall_spline], name="wall")
    gmsh.model.addPhysicalGroup(1, [inlet_line],  name="inlet")
    gmsh.model.addPhysicalGroup(1, [exit_line],   name="outlet")
    gmsh.model.addPhysicalGroup(1, [axis_line],   name="axis")
    gmsh.model.addPhysicalGroup(2, [surface],     name="fluid")

    # ── Mesh and export ───────────────────────────────────────────────────
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")      # smooth internal nodes

    gmsh.write("nozzle.msh")                   # GMSH format (for GUI inspection)
    gmsh.write("nozzle.su2")                   # SU2 format (for solver)

    print("Mesh written to nozzle.msh and nozzle.su2")
    print(f"Nodes:    {gmsh.model.mesh.getNodes()[0].shape[0] // 3}")

    gmsh.finalize()

if __name__ == "__main__":
    main()
```

**Notes on the GMSH script:**
- The `BoundaryLayer` field in GMSH requires GMSH version ≥ 4.8. Install the latest stable release.
- The physical group names (`wall`, `inlet`, `outlet`, `axis`) **must match** the `MARKER_*` names
  in the SU2 config file exactly — this is a common source of setup errors.
- For the Euler (inviscid) case, remove the BoundaryLayer field and use a simpler uniform triangle
  mesh — no boundary layer resolution is needed since the wall is a slip boundary.
- After generating the mesh, open `nozzle.msh` in the GMSH GUI to visually inspect the cell quality
  near the throat before running SU2.

---

## Benchmarking Against Known Nozzles

Before validating a custom design, validate the **MOC code itself** against published data. If the
MOC code cannot reproduce known results, there is no point in proceeding to CFD validation of a
custom nozzle.

### Rao's 1958 Nozzle Data

G.V.R. Rao's 1958 paper *"Exhaust Nozzle Contour for Optimum Thrust"* (Jet Propulsion, Vol. 28)
contains wall contour data tables for optimum nozzles at various area ratios and M_exit values. The
Rao nozzles are different from the uniform-flow nozzles produced by this MOC code (Rao optimizes
for maximum thrust in a fixed nozzle length rather than uniform exit flow), but they are well-known
reference data.

**What to compare:**
- For a given A_e/A_t and M_exit, does the MOC code produce a contour that gives approximately the
  same wall length as Rao's optimum nozzle?
- A uniform-flow (maximum uniformity) nozzle will generally be longer than a Rao optimum nozzle
  for the same area ratio — verify this is the case.

### NASA SP-8120 (1976)

*"Liquid Rocket Engine Nozzles"*, NASA SP-8120 (1976), provides performance data for various nozzle
types including:
- Theoretical C_F and Isp as a function of P_c/P_e and area ratio
- Correction factors for conical nozzles
- Boundary layer loss estimates

Compare your MOC-predicted C_F against the SP-8120 tables for the same γ, A_e/A_t, and P_c/P_e.
Discrepancies larger than 2% warrant investigation.

### Simple Self-Validation Test

A basic check that requires no external data:

1. Design a **conical nozzle** with half-angle θ = 15° and area ratio A_e/A_t = 8 using the
   simple geometric formula
2. Use the MOC code to design a **uniform-exit-flow nozzle** for the same A_e/A_t and M_exit
3. Run both contours through the **Euler CFD** simulation
4. Compare exit Mach uniformity: the MOC contour should show < 1% variation across the exit plane,
   while the conical nozzle should show 5–15% variation (the divergence loss)
5. This comparison validates that the MOC code does what it is supposed to do — eliminate flow
   non-uniformity — even if you cannot compare the absolute numbers to published data

### von Kármán Institute (VKI) Nozzle Database

The Von Kármán Institute for Fluid Dynamics (Brussels) has published experimental and numerical data
for several supersonic nozzle configurations. Their lecture series notes on *"Design and Testing of
High Performance Nozzles"* (various years) contain detailed wall pressure and Mach number profiles
useful for benchmarking both MOC codes and CFD setups.

---

## Summary

- **CFD validates** that the MOC-designed contour actually produces shock-free, uniform exit flow in
  a real viscous simulation — it is the final engineering check before hardware is committed.
- **Recommended tool:** SU2 (free, axisymmetric Euler + RANS, well-documented, actively maintained);
  OpenFOAM is a strong alternative with a steeper learning curve.
- **Two-step approach:** run **Euler first** to validate the contour shape (check for shocks),
  then add **RANS** to quantify viscous losses and compare with the V4 boundary layer estimate.
- **Workflow:** export contour from Rust as CSV → mesh in GMSH with BL refinement → solve in SU2
  → visualize and extract metrics in ParaView.
- **Key comparison metrics:** exit Mach uniformity (< 2% std dev), wall pressure profile, absence
  of shocks in the divergent section, thrust coefficient C_F within 2% of isentropic prediction.
- **Benchmark first:** validate the MOC code against Rao (1958) or the self-validation conical vs.
  uniform-flow comparison before running custom designs.
- **V5 is not a single implementation** but an ongoing validation practice: re-run the CFD
  comparison every time the MOC solver, gas model (V2/V3), or boundary layer correction (V4) is
  updated to confirm that improvements are real and regressions are caught early.
