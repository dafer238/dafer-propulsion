# V1 — Axisymmetric Method of Characteristics

This document covers the first and highest-priority enhancement to the MOC
nozzle solver: replacing the 2D planar geometry with the physically correct
axisymmetric geometry. It derives the modified governing equations from first
principles, explains each modification in physical terms, and gives the exact
Rust code changes required.

---

## Why Axisymmetric?

The baseline solver (V0) is a 2D **planar** solver. "Planar" means the flow is
identical in every plane parallel to the x-y plane — the nozzle is effectively
infinite in the z-direction, like a rectangular slit nozzle. The governing
equations are derived in 2D Cartesian coordinates (x, y), and the results are
a 2D streamline pattern in the x-y plane.

Real rocket nozzles are **axisymmetric**: they have circular cross-sections,
and the flow has no variation in the azimuthal (φ) direction. The governing
equations are derived in cylindrical coordinates (x, r), where r is the radial
distance from the centerline.

This is not merely a geometric labeling difference. The physics changes
fundamentally:

### Area Growth Rate

| Geometry | Cross-sectional Area | Growth with y or r |
|----------|---------------------|--------------------|
| Planar | A ∝ y (width × unit depth) | Linear |
| Axisymmetric | A = π r² | Quadratic |

For the same local half-angle θ of the wall, the divergence of the flow is
twice as strong in the axisymmetric case (the streamlines diverge both upward
and inward simultaneously). This stronger divergence:

- Accelerates the flow more rapidly for a given wall angle
- Means the flow needs a steeper wall to reach the same exit Mach number
- Changes the Prandtl-Meyer turning at every interior point

### The Radial Divergence Source Term

When you write Newton's second law and the continuity equation in cylindrical
coordinates for axisymmetric flow and manipulate them into the form of
characteristic equations, a new term appears that is not present in the planar
case. This term is proportional to v/r, where v is the radial velocity component
and r is the distance from the axis.

Physically, this term represents the fact that even a purely radial flow (θ > 0)
causes the cross-sectional area to grow as the flow moves radially outward,
which is equivalent to a source of mass in the planar interpretation. The
planar equations do not account for this, so they underpredict the turning that
any given nozzle section provides to the flow.

### Visible Effect on the Nozzle Contour

For area ratios A_e/A* ≤ 3, the planar and axisymmetric contours are similar
because the nozzle is short and the radial coordinates remain modest. For
A_e/A* > 5:

- The axisymmetric wall is **steeper** near the throat (more aggressive
  initial expansion)
- The axisymmetric nozzle is **physically shorter** for the same A_e/A*
  (the stronger divergence reaches M_exit sooner in x)
- The exit wall angle is **larger** in the axisymmetric case

For A_e/A* = 10 (the example in this project), overlaying the two contours on a
single plot will show a clearly visible difference in wall shape. The
axisymmetric contour is not a small correction — it is the correct contour.

**No real rocket nozzle is planar.** V1 is not an optional refinement; it is a
correction to a fundamental geometric error.

---

## The Governing Equation for Axisymmetric Flow

### Starting Point: The Velocity Potential Equation

For steady, irrotational, isentropic flow, the velocity can be written as the
gradient of a potential: **q** = ∇φ, giving u = ∂φ/∂x and v = ∂φ/∂r. The
governing equation for φ in cylindrical coordinates (x, r) for axisymmetric
supersonic flow is:

```
(1 - u²/a²) φ_xx - 2(uv/a²) φ_xr + (1 - v²/a²) φ_rr + (1/r)(1 - v²/a²) φ_r = 0
```

where subscripts denote partial derivatives and a is the local speed of sound.
This can be rewritten in terms of the velocity components directly. Using
u = q cos(θ), v = q sin(θ), M = q/a, and collecting terms:

```
(1 - u²/a²) φ_xx - 2(uv/a²) φ_xr + (1 - v²/a²) φ_rr - (v/r) = 0
```

The last term, `−v/r`, is the **axisymmetric source term**. It has no
counterpart in 2D planar flow. Note:

- When v = 0 (parallel flow, θ = 0), the source term vanishes. The flow along
  the nozzle centerline (where θ = 0) sees no axisymmetric correction.
- When r → ∞, the source term → 0. Far from the axis, the local geometry
  becomes planar again.
- The source term is largest where v is large (steep wall angle) and r is small
  (near the throat and near the axis). This is exactly where most of the
  turning takes place, which is why the correction matters.
- The source term is negative (it subtracts from the LHS), meaning it acts to
  **reduce** the effective angle that each interior point "sees", which requires
  the wall to turn more steeply to achieve the same expansion.

### The Characteristic Directions

The characteristic directions of the axisymmetric velocity potential equation
are the same as in the planar case:

```
dy/dx = tan(θ ± μ)
```

where μ = arcsin(1/M) is the Mach angle. **The characteristic slopes are
unchanged by the axisymmetric term.** This is important: you can reuse the
position calculation from V0 without modification. The only change is in the
compatibility equations that say what happens to the flow properties along those
characteristics.

---

## Modified Compatibility Equations

### Planar Case (V0 Recap)

In 2D planar MOC, the Riemann invariants K± are **exactly** constant along
their respective characteristics:

```
K⁺ = θ + ν(M)   is constant along C⁻   [slope = tan(θ - μ)]
K⁻ = θ - ν(M)   is constant along C⁺   [slope = tan(θ + μ)]
```

where ν(M) is the Prandtl-Meyer function. The beauty of planar MOC is that
you can find the state at any point simply by intersecting two known invariants.

### Axisymmetric Case (V1)

In axisymmetric MOC, the Riemann invariants are **not** constant along the
characteristics. The axisymmetric source term continuously modifies them as
you travel along each characteristic. The compatibility equations become
ordinary differential equations rather than algebraic statements:

**Along C⁺** (characteristic with slope = tan(θ + μ), which carries information
about K⁻):

```
d(θ - ν) / dx = + [sin(μ) · sin(θ)] / [r · cos(θ + μ)]
```

**Along C⁻** (characteristic with slope = tan(θ - μ), which carries information
about K⁺):

```
d(θ + ν) / dx = - [sin(μ) · sin(θ)] / [r · cos(θ - μ)]
```

where r is the **radial distance from the axis** — this is the `y` coordinate
in the existing code's x-y plane.

### Physical Interpretation of Each Term

- **sin(μ)**: scales with 1/M for large M (small Mach angles mean the source
  term diminishes at high Mach numbers in the exit region)
- **sin(θ)**: zero at θ = 0 (axis), maximum where the wall angle is steepest
- **r in the denominator**: large near the throat/axis (where corrections are
  largest), small far from the axis
- **cos(θ + μ) or cos(θ - μ)**: close to 1 for moderate angles, providing
  only a mild modulation
- **Positive for C⁺, negative for C⁻**: the source term acts to increase θ
  along C⁺ characteristics and to increase ν (decrease θ) along C⁻. The net
  effect is that the axisymmetric expansion is more efficient — the flow turns
  more for the same expansion, or equivalently, less wall angle is needed to
  achieve a given Mach number. However, the contour comes out steeper overall
  because the exit area ratio constraint forces more total turning.

### Key Observations for Implementation

1. **Axis singularity (r = 0)**: The source terms contain r in the denominator.
   At the axis, r = 0, but θ = 0 as well (the centerline is a streamline, so
   the flow must be parallel there). The product sin(θ)/r = sin(θ)/r → 0/0.
   By L'Hôpital's rule, this limit is finite and can be evaluated. In practice,
   because axis nodes in the code have θ = 0 exactly, the numerator is zero and
   the source term is zero. See the dedicated section below.

2. **Large-r limit**: For nodes far from the axis, the source term is small.
   For a nozzle with a large throat radius (r_t = 100 mm), points 200 mm from
   the axis see a source term that is ~50× smaller than points 4 mm from the
   axis. This is why planar and axisymmetric results converge for geometrically
   large nozzles with small divergence angles.

3. **Positive source on C⁺, negative on C⁻**: The source term always acts in
   the direction that increases θ (the flow direction angle) compared to the
   planar solution. This makes physical sense: the radial divergence forces the
   flow to turn more aggressively toward the axis.

---

## The Finite-Difference Approximation

To integrate the compatibility ODEs numerically, we replace the differential
along each characteristic with a finite difference over the interval from the
known upstream node to the unknown downstream point P.

### Setup for the Interior Point

Given lower node L (on C⁺) and upper node R (on C⁻), the interior point P is
found as follows.

**Step 1 — Position of P** (unchanged from V0):

```
x_P = [y_R - y_L + tan(θ_L + μ_L)·x_L - tan(θ_R - μ_R)·x_R]
       ÷ [tan(θ_L + μ_L) - tan(θ_R - μ_R)]

y_P = y_L + tan(θ_L + μ_L) · (x_P - x_L)
```

The position calculation is exactly the same as in V0. The characteristic
slopes do not change.

**Step 2 — Integrate the compatibility equations from L to P and R to P**:

The source term contributions (first-order in Δx):

```
ΔK⁻ from L to P  =  + [sin(μ_L) · sin(θ_L)] / [y_L · cos(θ_L + μ_L)] · (x_P - x_L)

ΔK⁺ from R to P  =  - [sin(μ_R) · sin(θ_R)] / [y_R · cos(θ_R - μ_R)] · (x_P - x_R)
```

Note that r = y in the solver's coordinate system (y is the distance from the
axis).

**Step 3 — Effective invariants arriving at P**:

```
K⁻_eff = (θ_L - ν_L) + ΔK⁻
K⁺_eff = (θ_R + ν_R) + ΔK⁺
```

**Step 4 — State at P**:

```
θ_P = (K⁺_eff + K⁻_eff) / 2
ν_P = (K⁺_eff - K⁻_eff) / 2
M_P = ν⁻¹(ν_P)          [inverse Prandtl-Meyer function]
```

This is the same algebraic structure as V0; only the effective invariants have
changed.

### The Predictor-Corrector Problem

The source term contributions ΔK⁻ and ΔK⁺ use x_P, y_P (via y_L), and the
states at L and R. The position x_P is already computed (Step 1 is unchanged).
However, a more accurate approximation would evaluate the source term using the
average state along the characteristic — i.e., using some estimate of the state
at P itself, which is what we're computing. This creates a chicken-and-egg
problem.

**First-order approximation (recommended for V1)**:

Use only the upstream node values (L or R) for the source term coefficients.
x_P from Step 1 is used as the Δx. This is first-order accurate and introduces
an error proportional to (Δx)² in the characteristic integrals. For a
well-resolved mesh (large n_chars), this error is small.

**Predictor-corrector (optional improvement)**:

1. Compute x_P, y_P from Step 1 (2D predictor position)
2. Compute ΔK⁻ and ΔK⁺ using the upstream nodes (first-order)
3. Solve for θ_P, ν_P, M_P (predicted state)
4. Re-evaluate the source terms using the average of (L, P) and (R, P)
5. Recompute the position x_P using average slopes
6. Re-solve for the final θ_P, ν_P, M_P

The corrector step improves accuracy from first-order to second-order in Δx.
For most nozzle design purposes with n_chars ≥ 30, the first-order approximation
is sufficient. For research-grade accuracy or very coarse meshes, use the
predictor-corrector.

### Axis Singularity in Finite Differences

When processing a characteristic that starts or passes near the axis (y_L < ε
or y_R < ε), the source term has y in the denominator. The safe numerical
treatment:

- If y_node < ε_r, set the source term for that node to 0.0
- The typical threshold is ε_r = 1×10⁻⁸ meters (or 1×10⁻⁶ × r_throat)
- Axis nodes always have θ = 0 exactly, so sin(θ) = 0 and the source term is
  analytically zero regardless of r; the threshold check is just a guard
  against floating-point division

---

## Rust Implementation

The implementation requires four changes to the codebase:

1. Add `axisymmetric: bool` to `NozzleConfig`
2. Add two source-term helper functions in `moc/characteristics.rs`
3. Modify `interior_point`, `axis_point`, and `wall_point` to add the
   source terms
4. Pass the flag through `MocSolver::design()`

### Step 1: Add `axisymmetric` Flag to `NozzleConfig`

In whatever file defines `NozzleConfig` (e.g., `src/config.rs` or `src/lib.rs`):

```rust
pub struct NozzleConfig {
    pub gamma:         f64,
    pub ae_at:         f64,
    pub n_chars:       usize,
    pub throat_radius: f64,
    pub axisymmetric:  bool,   // NEW: true = axisymmetric, false = 2D planar
}

impl Default for NozzleConfig {
    fn default() -> Self {
        Self {
            gamma:         1.4,
            ae_at:         10.0,
            n_chars:       30,
            throat_radius: 0.025,
            axisymmetric:  true,   // Default to physically correct geometry
        }
    }
}
```

Setting the default to `true` is appropriate because all real nozzles are
axisymmetric. Tests that want to reproduce the V0 planar behavior can
explicitly set `axisymmetric: false`.

### Step 2: Source Term Helpers in `moc/characteristics.rs`

Add these two functions near the top of the characteristics module, before the
node-computation functions:

```rust
/// Source term contribution to K⁻ (the C⁺ characteristic invariant)
/// from a node `l` to an estimated downstream position `x_p`.
///
/// Returns 0.0 when:
///   - axisymmetric is false (planar mode)
///   - the node is at or near the axis (y < 1e-8), where sin(θ) → 0
///   - the denominator is degenerate
///
/// The sign convention matches Anderson (1990) eq. 11.28:
///   dK⁻ = + sin(μ)·sin(θ) / [r · cos(θ + μ)] · dx
pub fn source_km(l: &Node, x_p: f64, axisymmetric: bool) -> f64 {
    if !axisymmetric { return 0.0; }
    if l.y < 1.0e-8 { return 0.0; }   // axis or near-axis: sin(θ) ≈ 0

    let mu    = (1.0 / l.state.m).asin();
    let theta = l.state.theta;
    let denom = l.y * (theta + mu).cos();

    if denom.abs() < 1.0e-12 { return 0.0; }

    mu.sin() * theta.sin() / denom * (x_p - l.x)
}

/// Source term contribution to K⁺ (the C⁻ characteristic invariant)
/// from a node `r` to an estimated downstream position `x_p`.
///
/// Returns 0.0 when:
///   - axisymmetric is false (planar mode)
///   - the node is at or near the axis
///   - the denominator is degenerate
///
/// The sign convention:
///   dK⁺ = - sin(μ)·sin(θ) / [r · cos(θ - μ)] · dx
pub fn source_kp(r: &Node, x_p: f64, axisymmetric: bool) -> f64 {
    if !axisymmetric { return 0.0; }
    if r.y < 1.0e-8 { return 0.0; }

    let mu    = (1.0 / r.state.m).asin();
    let theta = r.state.theta;
    let denom = r.y * (theta - mu).cos();

    if denom.abs() < 1.0e-12 { return 0.0; }

    -mu.sin() * theta.sin() / denom * (x_p - r.x)
}
```

#### Why Two Separate Functions?

`source_km` is applied to the **lower** (C⁺) node when computing K⁻ at P.
`source_kp` is applied to the **upper** (C⁻) node when computing K⁺ at P.
They differ only in the sign and in whether they use `cos(θ + μ)` or
`cos(θ - μ)`. Keeping them separate avoids sign confusion at the call site.

### Step 3: Modify `interior_point`

The updated `interior_point` function with the axisymmetric flag. Compare this
carefully with the V0 version — only the K± computation changes:

```rust
pub fn interior_point<G: GasModel>(
    lower: &Node,
    upper: &Node,
    gas: &G,
    axisymmetric: bool,
) -> Node {
    let mu_l = gas.mach_angle(lower.state.m);
    let mu_u = gas.mach_angle(upper.state.m);

    // --- Position of P (unchanged from V0) ---
    let s_plus  = (lower.state.theta + mu_l).tan();   // C⁺ slope from L
    let s_minus = (upper.state.theta - mu_u).tan();   // C⁻ slope from R
    let denom   = s_plus - s_minus;

    let x_p = if denom.abs() < 1.0e-12 {
        (lower.x + upper.x) / 2.0
    } else {
        (upper.y - lower.y + s_plus * lower.x - s_minus * upper.x) / denom
    };
    let y_p = lower.y + s_plus * (x_p - lower.x);

    // --- Invariants with axisymmetric source terms ---
    // K⁻ arrives from L along C⁺; pick up the source term ΔK⁻
    let k_minus_base = lower.state.theta - lower.state.nu;
    let k_minus = k_minus_base + source_km(lower, x_p, axisymmetric);

    // K⁺ arrives from R along C⁻; pick up the source term ΔK⁺
    let k_plus_base = upper.state.theta + upper.state.nu;
    let k_plus = k_plus_base + source_kp(upper, x_p, axisymmetric);

    // --- Solve for state at P ---
    let theta_p = 0.5 * (k_plus + k_minus);
    let nu_p    = 0.5 * (k_plus - k_minus);
    let m_p     = gas.mach_from_nu(nu_p);

    Node {
        x: x_p,
        y: y_p,
        state: FlowState { theta: theta_p, nu: nu_p, m: m_p },
    }
}
```

The structure is identical to V0; the only new lines are the `source_km` and
`source_kp` calls on the invariants.

### Step 4: Modify `axis_point`

The axis point is the special case where the C⁺ characteristic reaches the
nozzle centerline (y = 0). At the axis, by symmetry, θ = 0. The C⁺ arrives
from lower node L (which is not on the axis) and the axis condition provides
the other constraint.

The axis constraint in the planar case is: θ_P = 0, and K⁺ = K⁻ (symmetric).
In the axisymmetric case, the same constraint applies at the axis (θ = 0 is
enforced by symmetry), but the arriving K⁻ from L picks up the source term:

```rust
pub fn axis_point<G: GasModel>(
    lower: &Node,
    gas: &G,
    axisymmetric: bool,
) -> Node {
    let mu_l = gas.mach_angle(lower.state.m);

    // The C⁺ characteristic from L hits the axis at y = 0
    // Slope: tan(θ_L + μ_L); extrapolate to y = 0
    let slope = (lower.state.theta + mu_l).tan();
    let x_p   = lower.x - lower.y / slope;   // y_P = 0 → Δx = -y_L / slope
    let y_p   = 0.0;

    // K⁻ from L along C⁺, with source term (using x_P as the destination)
    // Note: near the axis, lower.y may be small but is not zero (L is
    // interior), so source_km is evaluated with lower.y in the denominator.
    let k_minus_base = lower.state.theta - lower.state.nu;
    let k_minus = k_minus_base + source_km(lower, x_p, axisymmetric);

    // Axis condition: θ_P = 0, so K⁺ = -K⁻ at the axis
    // (From θ = (K⁺+K⁻)/2 = 0 → K⁺ = -K⁻)
    // Actually: K⁻ = θ - ν and K⁺ = θ + ν; at axis θ = 0:
    //   K⁻ = -ν_P  and  K⁺ = +ν_P
    // So ν_P = -K⁻ (if K⁻ is negative, ν_P is positive as expected)
    let nu_p = -k_minus;
    let m_p  = gas.mach_from_nu(nu_p);

    Node {
        x: x_p,
        y: y_p,
        state: FlowState { theta: 0.0, nu: nu_p, m: m_p },
    }
}
```

**Why is the source term non-zero here?** The source term for `source_km(lower,
x_p, axisymmetric)` uses the position of node L (not P). L is an interior node
with y_L > 0 and θ_L > 0 in general, so the source term is non-zero and
contributes a modification to the K⁻ invariant arriving at the axis. The axis
itself (P) has θ = 0 and r = 0, so any source term evaluated at P would be
zero, but we are evaluating the accumulated source term along the characteristic
from L to P, which passes through non-zero radii.

### Step 5: Modify `wall_point`

The wall point requires finding where a C⁻ characteristic from interior node J
intersects the nozzle wall. The wall imposes a boundary condition: the flow at
the wall must be parallel to the wall (θ_wall = wall_angle at that x).

```rust
pub fn wall_point<G: GasModel>(
    j_node: &Node,         // Interior node the C⁻ comes from
    prev_wall: &Node,      // Previous wall node (provides wall direction)
    wall_angle: f64,       // Target wall flow angle at the new point
    gas: &G,
    axisymmetric: bool,
) -> Node {
    let mu_j = gas.mach_angle(j_node.state.m);

    // C⁻ slope from J (the characteristic heading toward the wall)
    let slope_char = (j_node.state.theta - mu_j).tan();

    // Wall slope: average direction between prev_wall and the new wall point
    // (In practice, using prev_wall direction is a good first approximation)
    let slope_wall = (prev_wall.state.theta).tan();

    // Intersect C⁻ from J with the wall
    let denom = slope_char - slope_wall;
    let x_p = if denom.abs() < 1.0e-12 {
        j_node.x + 0.01 * j_node.x   // fallback: small step forward
    } else {
        (prev_wall.y - j_node.y + slope_char * j_node.x
            - slope_wall * prev_wall.x) / denom
    };
    let y_p = j_node.y + slope_char * (x_p - j_node.x);

    // K⁺ arriving from J along C⁻, with source term
    let k_plus_base = j_node.state.theta + j_node.state.nu;
    let k_plus = k_plus_base + source_kp(j_node, x_p, axisymmetric);

    // Wall boundary condition: θ_P = wall_angle
    // From K⁺ = θ + ν: ν_P = K⁺ - θ_P
    let nu_p = k_plus - wall_angle;
    let m_p  = gas.mach_from_nu(nu_p);

    Node {
        x: x_p,
        y: y_p,
        state: FlowState { theta: wall_angle, nu: nu_p, m: m_p },
    }
}
```

The wall point uses only `source_kp` (from J along C⁻). There is no C⁺
characteristic at the wall point — the wall itself provides the second
constraint (θ = wall_angle).

### Step 6: Pass the Flag Through `MocSolver::design()`

Wherever `interior_point`, `axis_point`, and `wall_point` are called in the
main design loop, pass `config.axisymmetric`:

```rust
// Interior point: replace
let p = interior_point(&lower, &upper, &gas);
// With:
let p = interior_point(&lower, &upper, &gas, config.axisymmetric);

// Axis point: replace
let p = axis_point(&lower, &gas);
// With:
let p = axis_point(&lower, &gas, config.axisymmetric);

// Wall point: replace
let p = wall_point(&j_node, &prev_wall, wall_angle, &gas);
// With:
let p = wall_point(&j_node, &prev_wall, wall_angle, &gas, config.axisymmetric);
```

This is a mechanical change with no logic other than threading the flag through.

---

## What Changes in the Nozzle Contour

Running the solver with `axisymmetric: true` vs. `axisymmetric: false` (all
other parameters equal) will produce two different wall contours. The
differences are predictable and can be used as a sanity check:

### Expected Differences (Axisymmetric vs. Planar, same A_e/A*, n_chars, γ)

| Property | Planar (V0) | Axisymmetric (V1) | Direction |
|----------|-------------|-------------------|-----------|
| Initial wall angle (near throat) | Smaller | Larger | Axi steeper |
| Nozzle length (x_exit - x_throat) | Longer | Shorter | Axi shorter |
| Exit wall angle | Smaller | Larger | Axi steeper |
| y_exit / r_throat | Matches √(A_e/A*) for planar | Matches √(A_e/A*) for axi | Both correct for their geometry |

### Validation Check

After running V1, validate the exit area ratio by computing:

```rust
let computed_ae_at = (wall.last().y / config.throat_radius).powi(2);
let error = (computed_ae_at - config.ae_at).abs() / config.ae_at;
assert!(error < 0.01, "Exit area ratio error > 1%: got {}, expected {}",
        computed_ae_at, config.ae_at);
```

This check ensures that the characteristic mesh has expanded to exactly the
right area ratio. If this fails, there is a bug in either the source term signs
or the K± reconstruction.

### Cross-Validation Against Published Data

For γ = 1.4 and A_e/A* = 10.0, compare the computed exit Mach number against
the isentropic table:

```
M_exit = ν⁻¹(ν_exit)   where ν_exit satisfies A_e/A*(M_exit, γ=1.4) = 10.0
```

The isentropic M_exit for A_e/A* = 10.0 at γ = 1.4 is approximately 3.96.
All wall nodes near the exit should have θ ≈ 0 (parallel flow) and M ≈ 3.96.
If the exit Mach numbers from the characteristic network are far from this
value, the source terms have the wrong sign or magnitude.

---

## The Axis Singularity: L'Hôpital's Rule

This section provides the mathematical justification for why the source term
can be set to zero at axis nodes.

### The Limit

At a node exactly on the axis: y = 0 and θ = 0. The source term in the
compatibility equations is:

```
sin(μ) · sin(θ) / [y · cos(θ + μ)]
```

At y = 0, θ = 0:

```
= sin(μ) · sin(0) / [0 · cos(0 + μ)]
= sin(μ) · 0 / [0 · cos(μ)]
= 0 / 0   (indeterminate)
```

### Applying L'Hôpital's Rule

Consider the limit as a streamline approaches the axis. As r → 0, the
streamline angle θ → 0 as well (the flow approaches the centerline
asymptotically). The ratio:

```
lim_{r→0, θ→0} sin(θ) / r
```

equals `dθ/dr` evaluated at r = 0, which is a finite quantity related to the
local flow curvature. Specifically, for a smooth nozzle, the centerline is a
line of symmetry and θ changes continuously from zero on the axis. So
`sin(θ)/r → θ/r → (dθ/dr)|_{r=0}` (using small-angle approximation for θ).

### Numerical Consequence

In the code, axis nodes are computed by `axis_point`, which enforces θ = 0.
When this node is subsequently used as an upstream node in `source_km` or
`source_kp`, its θ = 0 makes sin(θ) = 0, and the source term evaluates to 0.0
regardless of the y value. The guard `if l.y < 1e-8 { return 0.0; }` is purely
defensive — the result would be 0.0 even without it, as long as θ = 0.

The only numerical risk arises from finite precision: if a node is computed to
have y = 1e-10 (effectively on the axis) but θ = 1e-7 (not exactly zero due to
rounding), the source term `sin(θ)/y ≈ 1e-7 / 1e-10 = 1000` would be
spuriously large. The guard threshold eliminates this edge case.

---

## Summary

The axisymmetric correction is the single most important enhancement to the
baseline solver. Key points:

- The axisymmetric source term `sin(μ)·sin(θ) / [r·cos(θ±μ)] · Δx` must be
  added to both compatibility equations at every interior, axis, and wall node
- The source term is **positive for K⁻ (C⁺ characteristics)** and
  **negative for K⁺ (C⁻ characteristics)**, meaning it always acts to produce
  greater effective turning compared to the planar case
- The **position computation is unchanged** — characteristic slopes are the
  same in planar and axisymmetric flow; only the invariants change
- Use the upstream node's state and the computed (2D-predictor) position for
  evaluating the source term — first-order approximation, sufficient for
  n_chars ≥ 20
- The **axis singularity is handled trivially**: axis nodes have θ = 0 exactly,
  so sin(θ) = 0 and the source term is zero; a numerical guard prevents
  division by near-zero r from amplifying floating-point noise
- Implementation touches exactly **four places**: `NozzleConfig`, `source_km`,
  `source_kp`, and the three node functions (`interior_point`, `axis_point`,
  `wall_point`)
- **Every real rocket nozzle is axisymmetric** — this version is not optional
  if the goal is a nozzle contour that can be machined and tested

After implementing V1, validate against published Rao nozzle data before
proceeding to V2/V3.

---

*Previous: [07 — Enhancement Roadmap](07_enhancements_roadmap.md)*  
*Next: [09 — V2 Variable Specific Heats](09_v2_variable_gamma.md)*
