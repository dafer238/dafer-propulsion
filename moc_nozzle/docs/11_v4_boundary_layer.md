# V4 — Boundary-Layer Displacement Correction

## What a Boundary Layer Is

The MOC solver assumes **inviscid** flow — no friction, no viscosity. In reality, viscosity causes
the flow to stick to the nozzle wall (the **no-slip condition**). A thin region near the wall, called
the **boundary layer**, transitions from zero velocity at the wall to the full inviscid free-stream
velocity just outside it.

The flow profile across the boundary layer looks approximately like this (r is distance from wall):

```
u/u_e
 1.0  |                     ●●●●●●●●●●●●  ← free stream (inviscid core)
      |                ●●
      |             ●●
      |           ●●
 0.5  |         ●
      |        ●
      |       ●
      |      ●
 0.0  |●●●●●  ← no-slip wall
       ─────────────────────────────────► r (distance from wall)
              δ (boundary layer thickness)
```

The key consequence is that the boundary layer **displaces** the effective inviscid core away from
the wall. The **displacement thickness δ\*(x)** is defined as the thickness you would need to move
the wall inward to make the inviscid flow carry exactly the same mass flow as the real viscous flow.
Formally:

```
         ∞
δ*(x) = ∫  [ 1 - ρ(r)·u(r) / (ρ_e·u_e) ] dr
         0
```

For the nozzle, this has several practical consequences:

- At the throat, δ\* reduces the effective throat area: `A*_eff = π(r_t - δ*_t)²`
- This smaller effective throat area means the nozzle passes **less mass flow** than the inviscid
  prediction
- The exit area ratio `A_e / A*_eff` is **larger** than `A_e / A*_inviscid`
- The exit Mach number is slightly **lower** than the inviscid prediction
- For **small nozzles** (r_t < 10 mm, such as laboratory thrusters), this correction is critical —
  δ\*/r_t can be several percent and significantly shifts the achieved Mach number
- For **large nozzles** (r_t > 50 mm, typical orbital rockets), the correction is < 1% on Isp, but
  still affects the precise contour shape used in high-fidelity design

The boundary layer also generates a **displacement body**: the effective wall seen by the inviscid
core is shifted inward by δ\*(x) at every axial station. The corrected wall contour is therefore:

```
y_corrected(x) = y_MOC(x) - δ*(x)
```

If you want the nozzle to achieve exactly M_exit at the geometric wall, you must **increase** the
geometric wall radius by δ\*(x) everywhere — i.e., the physical hardware wall sits at
`y_hardware(x) = y_MOC(x) + δ*(x)`. The two perspectives (corrected contour vs. enlarged hardware)
are equivalent.

---

## Physics of the Compressible Boundary Layer

### Laminar Boundary Layer

For a laminar boundary layer on a flat plate in compressible flow (Blasius solution, corrected for
compressibility):

```
δ*(x) ≈ (1.72 / √Re_x) · x · I_lam
```

where:
- `Re_x = ρ_e · u_e · x / μ_e` is the local Reynolds number based on distance x from the
  throat
- `I_lam` is a compressibility integral that approaches 1.0 at low Mach and increases with M;
  for M < 3 it stays below ~1.3

The local edge conditions (subscript e) come directly from the MOC wall node states via isentropic
relations:

```
T_e(x) = T_0 / [ 1 + (γ-1)/2 · M(x)² ]
P_e(x) = P_0 · [ T_e / T_0 ]^(γ/(γ-1))
ρ_e(x) = P_e / (R_specific · T_e)
u_e(x) = M(x) · √(γ · R_specific · T_e)
```

### Turbulent Boundary Layer

Most rocket nozzles operate at Reynolds numbers above 10⁵ at the throat, which means the boundary
layer transitions to turbulent early in the diverging section. Turbulent boundary layers are thicker
and grow faster than laminar ones.

The **Eckert reference temperature method** is the standard engineering approach for compressible
turbulent boundary layers:

**Step 1 — Compute reference temperature T\*:**

```
T*(x) = 0.5·(T_w + T_e) + 0.22·r_f·[(γ-1)/2]·M_e²·T_e
```

where:
- `T_w` = wall temperature (typically 500–1000 K for regeneratively cooled nozzles, up to 3000 K
  for uncooled ablative walls)
- `T_e` = local free-stream temperature from MOC isentropic relations
- `r_f` = recovery factor ≈ Pr^(1/3) ≈ 0.89 for turbulent boundary layers
  (≈ √Pr ≈ 0.84 for laminar)
- `M_e` = local edge Mach number

The recovery factor r_f accounts for the fact that the adiabatic wall temperature is slightly higher
than T_e due to viscous heating (the wall "recovers" a fraction of the kinetic energy as heat).

**Step 2 — Compute reference state properties:**

```
ρ*(x) = P_e / (R_specific · T*)      ← same pressure, reference temperature
μ*(x) = μ_sutherland(T*)              ← viscosity at reference temperature
Re*_x = ρ*(x) · u_e(x) · x / μ*(x)  ← reference Reynolds number
```

**Step 3 — Turbulent skin friction at reference conditions:**

Using the Prandtl 1/7 power-law turbulent profile:

```
C_f*(x) = 0.0592 / Re*_x^(1/5)
```

This is valid for Re_x in the range 10⁵ to 10⁷. For higher Reynolds numbers, use the Schlichting
formula: `C_f* = 0.455 / [ln(0.06 · Re*_x)]²`

**Step 4 — Momentum and displacement thickness:**

From the von Kármán integral momentum equation (simplified for zero pressure gradient as an
approximation):

```
θ(x) = (C_f*(x) / 2) · x         ← momentum thickness [m]
δ*(x) = H* · θ(x)                 ← displacement thickness [m]
```

The turbulent shape factor H\* ≈ 1.4 for mild compressibility (M < 3). At higher Mach numbers,
H\* rises because the density near the wall drops (hot gas due to aerodynamic heating is less dense,
so you need more wall offset to recover the same mass defect).

**Simplified formula (good to within ~20% for engineering use):**

For turbulent BL on the nozzle wall at position x from the throat:

```
δ*(x) / r_t ≈ 0.664 · (ν_ref / (u_ref · r_t))^(1/5) · (x/r_t)^(4/5)
            = A · (x/r_t)^(4/5)
```

where A is a dimensional prefactor depending on gas properties and throat conditions.

For typical rocket nozzle conditions (P_t ~ 2 MPa, T_t ~ 2000 K, r_t = 25 mm):
- A ~ 0.002 to 0.005
- At x/r_t = 5 (exit of a typical nozzle): δ\*/r_t ~ 0.008 to 0.02

This means for a 25 mm throat radius nozzle, the boundary layer correction at the exit is on the
order of 0.2–0.5 mm — small but measurable, and relevant for Mach number accuracy at the 1% level.

---

## The Throat Correction

The boundary layer grows from zero at the nozzle entrance and is thickest at the exit. However, the
**throat correction** is disproportionately important because the throat sets the mass flow for the
entire nozzle.

At the throat (x = 0 in the convention where the throat is the origin), δ\* is small but non-zero
(the boundary layer has been growing through the convergent section upstream). In the divergent
section, the BL starts fresh from the throat and grows along x.

The **effective throat radius** is:

```
r_t_eff = r_t - δ*_throat
```

The **effective exit radius** is:

```
r_e_eff = r_e - δ*_exit
```

The effective area ratio becomes:

```
(A_e / A*)_eff = (r_e_eff / r_t_eff)²
               = (r_e - δ*_e)² / (r_t - δ*_t)²
```

This is slightly **larger** than the geometric ratio `(r_e / r_t)²` because δ\*_e > δ\*_t (the BL
is thicker at the exit than at the throat, and the denominator shrinks more than the numerator).

To design a nozzle that achieves a specific M_exit, you must account for this: the **geometric** area
ratio must be increased slightly above the inviscid value. Equivalently, the MOC contour computes the
inviscid wall, and you add δ\*(x) to the wall to get the physical hardware wall.

Numerically, for our example (P_t ~ 2 MPa, T_t ~ 2000 K, r_t = 25 mm, A_e/A_t = 8, design M = 3.3):
- δ\*_t ~ 0.05 mm → r_t_eff = 24.95 mm → throat area error of ~0.4%
- δ\*_e ~ 0.3 mm → r_e_eff = 55.7 mm (geometric r_e ≈ 56 mm)
- Actual area ratio = (55.7/24.95)² ≈ 4.98 instead of (56/25)² = 5.02
- This shifts the effective exit Mach by ΔM ≈ 0.015 (~0.5%)

For a 1 N thruster with r_t = 1.5 mm, the same calculation gives ΔM ≈ 0.05 (~1.5%), which is
significant for accurate performance prediction.

---

## Viscosity Models

To compute Reynolds numbers, you need dynamic viscosity μ(T). Two standard models:

### Sutherland's Law

This is the recommended model for engineering accuracy over a wide temperature range:

```
μ(T) = μ_ref · (T/T_ref)^(3/2) · (T_ref + S) / (T + S)
```

Standard constants for common gases:

| Gas        | μ_ref [Pa·s]          | T_ref [K] | S [K]  | Valid range     |
|------------|-----------------------|-----------|--------|-----------------|
| Air        | 1.716 × 10⁻⁵          | 273.15    | 110.4  | 170–1900 K      |
| H₂O vapor  | 1.12 × 10⁻⁵ (approx) | 373       | 673    | 373–1500 K      |
| CO₂        | 1.370 × 10⁻⁵          | 273.15    | 222.0  | 220–1800 K      |
| N₂         | 1.663 × 10⁻⁵          | 273.15    | 107.0  | 100–2000 K      |
| H₂         | 0.840 × 10⁻⁵          | 273.15    | 97.0   | 100–2000 K      |

For rocket exhaust (LOX/LH₂ mixture with M_mol ≈ 11 g/mol), use approximate values:
- μ_ref ≈ 1.0 × 10⁻⁵ Pa·s at T_ref = 1000 K, S ≈ 200 K

### Power Law

Simpler but less accurate, especially far from T_ref:

```
μ(T) = μ_ref · (T/T_ref)^n
```

where n ≈ 0.67 for diatomic gases (N₂, O₂, H₂) and n ≈ 0.69 for triatomic gases (CO₂, H₂O). The
power law is adequate when T stays within a factor of ~2 of T_ref.

### Mixture Viscosity for Rocket Exhaust

For frozen-flow rocket exhaust containing multiple species (H₂O, H₂, CO, CO₂, OH, etc.), use the
**Wilke mixing rule**:

```
μ_mix = Σ_i  x_i · μ_i / Σ_j  x_j · Φ_ij

Φ_ij = [1 + (μ_i/μ_j)^(1/2) · (M_j/M_i)^(1/4)]² / [8 · (1 + M_i/M_j)]^(1/2)
```

where x_i is the mole fraction of species i. This is complex to implement; for engineering
purposes, a **single effective viscosity** computed from CEA output at T_t is accurate to within
5–10% across the nozzle temperature range and is sufficient for BL calculations.

---

## Rust Implementation

Create a new file `src/geometry/boundary_layer.rs`:

```rust
use std::f64::consts::PI;
use crate::moc::node::Node;

/// Configuration for boundary layer calculation.
pub struct BoundaryLayerConfig {
    /// Wall temperature [K] — for regen cooling typically 500–1500 K,
    /// for uncooled ablative walls up to ~3000 K.
    pub t_wall: f64,

    /// Recovery factor (turbulent ≈ Pr^{1/3} ≈ 0.89, laminar ≈ √Pr ≈ 0.84).
    pub recovery: f64,

    /// Dynamic viscosity reference value [Pa·s] at t_ref (for Sutherland's law).
    pub mu_ref: f64,

    /// Reference temperature for Sutherland's law [K].
    pub t_ref: f64,

    /// Sutherland constant [K].
    pub sutherland_s: f64,

    /// Total (stagnation) pressure at chamber [Pa].
    pub p0: f64,

    /// Total (stagnation) temperature at chamber [K].
    pub t0: f64,

    /// Specific gas constant R_mix [J/(kg·K)] = R_universal / M_molar.
    pub r_specific: f64,

    /// Ratio of specific heats γ (at throat conditions).
    pub gamma_throat: f64,
}

impl BoundaryLayerConfig {
    /// Default configuration for LOX/LH₂ at the given chamber conditions.
    ///
    /// These values are approximate; for high-accuracy work, pull γ and R from
    /// your CEA output (V3).
    ///
    /// # Arguments
    /// * `p0` - Chamber stagnation pressure [Pa]
    /// * `t0` - Chamber stagnation temperature [K]
    pub fn lox_lh2(p0: f64, t0: f64) -> Self {
        Self {
            t_wall:       700.0,   // K, typical regen-cooled wall
            recovery:     0.89,    // turbulent
            mu_ref:       1.0e-5,  // Pa·s at t_ref = 1000 K
            t_ref:        1000.0,  // K
            sutherland_s: 200.0,   // K
            p0,
            t0,
            r_specific:   692.0,   // J/(kg·K), R_univ / 12 g/mol ≈ LOX/LH₂ mixture
            gamma_throat: 1.13,
        }
    }

    /// Configuration for a cold-gas N₂ thruster (lab testing / attitude control).
    pub fn nitrogen_cold_gas(p0: f64, t0: f64) -> Self {
        Self {
            t_wall:       t0,      // unheated wall ≈ gas temperature
            recovery:     0.89,
            mu_ref:       1.663e-5,
            t_ref:        273.15,
            sutherland_s: 107.0,
            p0,
            t0,
            r_specific:   296.8,   // J/(kg·K) for N₂
            gamma_throat: 1.40,
        }
    }
}

/// Dynamic viscosity [Pa·s] via Sutherland's law.
fn viscosity_sutherland(t: f64, cfg: &BoundaryLayerConfig) -> f64 {
    cfg.mu_ref
        * (t / cfg.t_ref).powf(1.5)
        * (cfg.t_ref + cfg.sutherland_s)
        / (t + cfg.sutherland_s)
}

/// Compute isentropic edge conditions at a given Mach number and stagnation state.
///
/// Returns (T_e, P_e, rho_e, u_e) where subscript e = "edge of boundary layer"
/// (i.e., the local free-stream condition from the inviscid MOC solution).
fn edge_conditions(
    mach: f64,
    gamma: f64,
    cfg: &BoundaryLayerConfig,
) -> (f64, f64, f64, f64) {
    let gm1 = gamma - 1.0;
    let factor = 1.0 + gm1 / 2.0 * mach * mach;

    let t_e   = cfg.t0 / factor;
    let p_e   = cfg.p0 * (t_e / cfg.t0).powf(gamma / gm1);
    let rho_e = p_e / (cfg.r_specific * t_e);
    let a_e   = (gamma * cfg.r_specific * t_e).sqrt();
    let u_e   = mach * a_e;

    (t_e, p_e, rho_e, u_e)
}

/// Compute displacement thickness δ*(x) at each wall node using the Eckert
/// reference temperature method for a turbulent boundary layer.
///
/// Assumes the boundary layer starts at the throat (x = 0). Nodes with x ≤ 0
/// return δ* = 0.
///
/// # Arguments
/// * `wall_nodes`  - Slice of wall nodes from the MOC solution (ordered x-ascending)
/// * `gamma`       - Effective γ (use the throat value)
/// * `cfg`         - Boundary layer configuration
///
/// # Returns
/// A `Vec<f64>` of δ* values [m], one per wall node, aligned with `wall_nodes`.
pub fn displacement_thickness(
    wall_nodes: &[Node],
    gamma:      f64,
    cfg:        &BoundaryLayerConfig,
) -> Vec<f64> {
    let gm1 = gamma - 1.0;

    wall_nodes.iter().map(|node| {
        let x = node.x;

        // No boundary layer at or upstream of the throat in this simplified model.
        if x <= 0.0 {
            return 0.0;
        }

        let m_e = node.state.m;

        // Isentropic edge conditions from MOC node state.
        let (t_e, p_e, _rho_e, u_e) = edge_conditions(m_e, gamma, cfg);

        // --- Eckert reference temperature ---
        // T* = 0.5*(T_w + T_e) + 0.22 * r_f * (γ-1)/2 * M_e² * T_e
        let t_star_raw = 0.5 * (cfg.t_wall + t_e)
            + 0.22 * cfg.recovery * gm1 / 2.0 * m_e * m_e * t_e;

        // Clamp to physical range: T* should be between T_e and T_wall in practice.
        let t_star = t_star_raw.clamp(t_e.min(cfg.t_wall), t_e.max(cfg.t_wall));

        // --- Reference state (pressure = local edge pressure, T = T*) ---
        let rho_star = p_e / (cfg.r_specific * t_star);
        let mu_star  = viscosity_sutherland(t_star, cfg);

        // Local Reynolds number at reference conditions.
        let re_x_star = rho_star * u_e * x / mu_star;

        if re_x_star < 1.0 {
            // Avoid singularity at the origin; return a small but nonzero value.
            return 0.0;
        }

        // --- Turbulent momentum thickness (Prandtl 1/7 power law) ---
        //   θ(x) = 0.0368 * x / Re*_x^(1/5)
        //
        // This integrates the skin friction C_f = 0.0592/Re_x^(1/5) over [0, x]
        // under the von Kármán momentum integral equation (zero pressure gradient
        // approximation).
        let theta = 0.0368 * x / re_x_star.powf(0.2);

        // --- Turbulent shape factor H* ---
        // H* = δ*/θ ≈ 1.4 for M < 3.
        // For higher Mach numbers a simple correction is H* ≈ 1.4 + 0.07*(M_e - 3)
        // for 3 < M_e < 6.
        let h_shape = if m_e <= 3.0 {
            1.4
        } else {
            1.4 + 0.07 * (m_e - 3.0)
        };

        theta * h_shape
    }).collect()
}

/// Apply the boundary layer correction: shift the wall contour inward by δ*(x).
///
/// Returns a `Vec<(f64, f64)>` of (x, y_corrected) pairs representing the effective
/// inviscid wall as seen by the inviscid core.
///
/// If you want the physical hardware wall instead, use `y_hardware = y_MOC + δ*(x)`.
///
/// # Arguments
/// * `wall_nodes`  - MOC wall nodes
/// * `delta_star`  - Displacement thicknesses from `displacement_thickness()`
pub fn corrected_wall(
    wall_nodes:  &[Node],
    delta_star:  &[f64],
) -> Vec<(f64, f64)> {
    wall_nodes
        .iter()
        .zip(delta_star.iter())
        .map(|(node, &ds)| {
            // Clamp to zero: the corrected wall cannot go below the axis.
            (node.x, (node.y - ds).max(0.0))
        })
        .collect()
}

/// Compute the effective area ratio correction due to the boundary layer.
///
/// Returns `(A_e_eff / A_t_eff)` — the area ratio the inviscid core actually "sees",
/// which is slightly larger than the geometric ratio.
///
/// # Arguments
/// * `r_throat`    - Geometric throat radius [m]
/// * `r_exit`      - Geometric exit radius [m]
/// * `ds_throat`   - δ* at the throat [m]
/// * `ds_exit`     - δ* at the exit [m]
pub fn effective_area_ratio(
    r_throat:  f64,
    r_exit:    f64,
    ds_throat: f64,
    ds_exit:   f64,
) -> f64 {
    let r_t_eff = r_throat - ds_throat;
    let r_e_eff = r_exit   - ds_exit;
    (r_e_eff / r_t_eff).powi(2)
}

/// Estimate the throat Reynolds number for transition-regime identification.
///
/// Re_throat ≈ ρ*_throat · u*_throat · r_t / μ*_throat
///
/// At the throat: M = 1, so all conditions follow directly from stagnation state.
pub fn throat_reynolds(r_throat: f64, gamma: f64, cfg: &BoundaryLayerConfig) -> f64 {
    let gm1 = gamma - 1.0;
    // Throat temperature (M=1): T* = T0 · 2/(γ+1)
    let t_throat  = cfg.t0 * 2.0 / (gamma + 1.0);
    let p_throat  = cfg.p0 * (t_throat / cfg.t0).powf(gamma / gm1);
    let rho_throat = p_throat / (cfg.r_specific * t_throat);
    let a_throat  = (gamma * cfg.r_specific * t_throat).sqrt(); // u = a at M=1
    let mu_throat = viscosity_sutherland(t_throat, cfg);

    rho_throat * a_throat * r_throat / mu_throat
}
```

### Module Registration

Add `pub mod boundary_layer;` to `src/geometry/mod.rs` alongside the existing geometry modules:

```rust
pub mod wall;
pub mod throat;
pub mod boundary_layer;   // ← add this line
```

### Integration in `main.rs`

```rust
use geometry::boundary_layer::{
    BoundaryLayerConfig,
    displacement_thickness,
    corrected_wall,
    effective_area_ratio,
    throat_reynolds,
};

// --- Chamber conditions ---
let p0     = 10.0e6;   // 10 MPa
let t0     = 3500.0;   // K
let gamma  = 1.13;
let r_throat = 0.025;  // 25 mm

// --- Boundary layer configuration ---
let bl_cfg = BoundaryLayerConfig::lox_lh2(p0, t0);

// --- Transition check ---
let re_t = throat_reynolds(r_throat, gamma, &bl_cfg);
println!("Throat Reynolds number: {:.2e}", re_t);
if re_t > 1.0e6 {
    println!("→ Turbulent BL assumed (Re > 10⁶)");
} else if re_t < 1.0e5 {
    println!("→ Laminar BL assumed (Re < 10⁵)");
} else {
    println!("→ Transitional regime — turbulent formula used conservatively");
}

// --- Compute displacement thickness at each wall node ---
let wall_nodes = nozzle.solver.wall_nodes();
let delta_star = displacement_thickness(wall_nodes, gamma, &bl_cfg);

// --- Apply correction to the wall contour ---
let corrected = corrected_wall(wall_nodes, &delta_star);

println!("=== Corrected wall contour (effective inviscid wall) ===");
println!("{:>10}  {:>12}  {:>12}  {:>12}", "x [m]", "y_MOC [m]", "δ* [m]", "y_eff [m]");
for ((x, y_eff), (node, ds)) in corrected.iter().zip(wall_nodes.iter().zip(delta_star.iter())) {
    println!("{:>10.4}  {:>12.6}  {:>12.6}  {:>12.6}", x, node.y, ds, y_eff);
}

// --- Effective area ratio ---
let r_exit    = wall_nodes.last().unwrap().y;
let ds_throat = delta_star[0];
let ds_exit   = *delta_star.last().unwrap();

let ae_at_geom = (r_exit / r_throat).powi(2);
let ae_at_eff  = effective_area_ratio(r_throat, r_exit, ds_throat, ds_exit);

println!("\nGeometric area ratio A_e/A_t = {:.4}", ae_at_geom);
println!("Effective area ratio (BL)    = {:.4}", ae_at_eff);
println!("Correction:                    {:.4}%", (ae_at_eff / ae_at_geom - 1.0) * 100.0);
```

---

## Laminar vs. Turbulent: When Does Each Apply?

| Condition | BL type | Typical Re_throat | δ\*/r_t at exit |
|-----------|---------|------------------|-----------------|
| Large orbital rocket (Vulcain 2, RS-25) | Fully turbulent | > 10⁷ | < 0.5% |
| Mid-scale bipropellant thruster (500 N) | Fully turbulent | ~10⁶ | 0.5–1% |
| Bipropellant thruster 100 N class | Transitional | 10⁵–10⁶ | 1–2% |
| Monopropellant thruster 1 N class | Mostly laminar | < 10⁵ | 3–8% |
| Cold-gas N₂ thruster | Laminar | < 10⁴ | > 10% |

**Rule of thumb:**
- Re_throat > 10⁶ → turbulent formulas throughout
- Re_throat < 10⁵ → laminar Blasius formula
- 10⁵ < Re_throat < 10⁶ → transitional; use turbulent formula as a conservative bound

**Quick estimate of Re_throat** (without a full flow-field calculation):

At the throat, M = 1 and the isentropic relations give:
```
T_throat = T0 · 2/(γ+1)
P_throat = P0 · [2/(γ+1)]^(γ/(γ-1))
ρ_throat = P_throat / (R · T_throat)
a_throat = √(γ · R · T_throat)            ← u_throat = a_throat at M=1
μ_throat = μ_sutherland(T_throat)

Re_throat = ρ_throat · a_throat · r_throat / μ_throat
```

For LOX/LH₂ at P₀ = 10 MPa, T₀ = 3500 K, r_t = 25 mm: Re_throat ≈ 4 × 10⁶ → turbulent.
For a 1 N cold-gas N₂ thruster at P₀ = 0.5 MPa, T₀ = 300 K, r_t = 1 mm: Re_throat ≈ 3 × 10⁴ →
laminar.

### Laminar Flat-Plate Formula (for completeness)

If you determine the BL is laminar, replace the turbulent momentum thickness formula with the
Blasius result:

```
θ_lam(x) = 0.664 · x / √Re_x_star       ← Blasius solution
H_lam     = 2.591                         ← shape factor for laminar BL
δ*_lam(x) = H_lam · θ_lam(x)
           = 1.721 · x / √Re_x_star
```

In Rust, this would replace the turbulent branch in `displacement_thickness()` with:

```rust
let theta   = 0.664  * x / re_x_star.sqrt();
let h_shape = 2.591;
theta * h_shape
```

---

## Summary

- The **boundary layer** forms because viscosity causes flow to stick to the wall (no-slip
  condition). The free-stream from the MOC is only an approximation valid outside this thin layer.
- The **displacement thickness δ\*(x)** quantifies how far the effective wall is shifted inward;
  the corrected wall is `y_eff(x) = y_MOC(x) − δ*(x)`.
- The most important effect is the **throat correction**: `r_t_eff = r_t − δ*_t` determines
  actual mass flow and, via the area ratio, the true exit Mach number.
- **Turbulent BL** dominates for Re_throat > 10⁶ (large rockets): use the Eckert reference
  temperature method.
- **Laminar BL** applies for Re_throat < 10⁵ (small cold-gas thrusters): use the Blasius flat-plate
  formula.
- **Viscosity**: Sutherland's law `μ(T) = μ_ref·(T/T_ref)^(3/2)·(T_ref+S)/(T+S)` is the
  recommended model.
- **Implementation**: compute δ\*(x) from wall node states (Mach number + isentropic relations),
  apply as a post-processing step on the MOC wall contour.
- **V4 is a post-processing step on V1+V2+V3**: no changes to the MOC kernel are needed; the
  correction is applied to the wall node output after the characteristic net is fully computed.
