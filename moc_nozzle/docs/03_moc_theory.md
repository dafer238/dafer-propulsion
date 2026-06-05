# Method of Characteristics — Theory

## Why MOC? The PDE Perspective

The steady, inviscid, irrotational supersonic flow in a nozzle is governed by a system of partial differential equations. The mathematical character of these PDEs depends on the local Mach number:

- When **M > 1** (supersonic): the system is **hyperbolic**. Information propagates along specific directions — the characteristic curves. The solution can be marched forward in space, from a known data line downstream toward the exit.
- When **M < 1** (subsonic): the system is **elliptic**. Information propagates in all directions simultaneously. The entire flow field must be solved globally — there is no marching direction, and MOC cannot be applied.

This is why MOC is the natural tool for the diverging (supersonic) section of a rocket nozzle, but cannot be used in the converging section. The throat marks the boundary between the two regimes.

For a hyperbolic system, **characteristics** are special curves in the solution domain along which the governing PDEs reduce to **ordinary differential equations** (the compatibility equations). This is a dramatic simplification: instead of solving coupled PDEs over a 2D domain, you integrate ODEs along a network of characteristic curves and advance the solution one point at a time.

The analogy with wave propagation makes this intuitive: in supersonic flow, a disturbance at a point can only influence the flow within the downstream Mach cone — a cone of half-angle μ = arcsin(1/M). Information cannot travel upstream. The characteristic lines are the boundaries of these Mach cones. Marching forward along them is exactly like tracing wavefronts in time-domain wave propagation, but unfolded into the spatial x-y plane.

This is in stark contrast with subsonic (elliptic) flow, where a disturbance at any point influences the entire flow field in all directions — just as a stone dropped in a still pond sends ripples everywhere, not just downstream.

---

## The Two Families of Characteristics

For 2D planar, steady, isentropic, supersonic flow, there are exactly two families of characteristic lines in the (x, y) plane.

### Notation

Throughout this project, the following notation is used consistently (all angles in **radians**):

| Symbol | Meaning |
|--------|---------|
| θ | Local flow angle — angle of the velocity vector from the x-axis |
| ν | Prandtl-Meyer angle — defined in `02_prandtl_meyer.md` |
| M | Local Mach number |
| μ | Mach angle = arcsin(1/M) |

### Characteristic Slopes and Invariants

The two characteristic families and their properties:

- **C⁺ characteristics** (left-running, "/" direction in the upper half-plane):
  - Slope: `dy/dx = tan(θ + μ)`
  - Carry the invariant: **K⁻ = θ − ν = constant**

- **C⁻ characteristics** (right-running, "\" direction in the upper half-plane):
  - Slope: `dy/dx = tan(θ − μ)`
  - Carry the invariant: **K⁺ = θ + ν = constant**

The naming convention (C⁺ carries K⁻, and vice versa) can be confusing at first. The superscript on K refers to the sign in the compatibility equation, not the characteristic family. It is a standard convention in MOC literature.

### Geometry of the Characteristic Network

```
Wall (y = r_wall)
   \          /
    \   C-   /  C+
     \      /
      \    /
       \  /
        \/  <- interior node P
        /\
       /  \
      /    \
   node R  node L
Axis (y = 0)
```

At interior node P:
- **Node L** (lower, closer to axis) lies on a C⁺ characteristic that reaches P — it provides **K⁻**.
- **Node R** (upper, closer to wall) lies on a C⁻ characteristic that reaches P — it provides **K⁺**.

Given K⁺ from node R and K⁻ from node L, the new state at P is fully determined:

```
θ_P = (K⁺ + K⁻) / 2
ν_P = (K⁺ − K⁻) / 2
M_P = inverse_prandtl_meyer(ν_P)
```

---

## Compatibility Equations

The compatibility equations are the ODEs that hold along each family of characteristics. They encode the conservation of the Riemann invariants:

**Along C⁺** (slope = tan(θ + μ)):

```
dθ − dν = 0   →   K⁻ = θ − ν = constant
```

**Along C⁻** (slope = tan(θ − μ)):

```
dθ + dν = 0   →   K⁺ = θ + ν = constant
```

This is what the code in `moc/characteristics.rs` implements:

```rust
pub fn invariants(s: FlowState) -> Invariants {
    Invariants {
        k_plus:  s.theta + s.nu,   // K⁺, constant along C⁻ characteristics
        k_minus: s.theta - s.nu,   // K⁻, constant along C⁺ characteristics
    }
}
```

This algebraic structure is correct. The bug is not here — it is that M is never recovered from ν after `from_invariants` computes the new θ and ν at a point. The new Mach number must be computed via `gas.inverse_prandtl_meyer(nu)`, which in turn requires `prandtl_meyer()` to return radians (not degrees, as it currently does).

---

## The Three Node Types

A MOC mesh is built by computing three types of nodes. Each type uses the same invariant algebra but with different geometric boundary conditions.

### Interior Point

An interior point P is computed from a left node L (on a C⁺ char) and a right node R (on a C⁻ char).

**Flow state at P:**

```
K⁻_L = θ_L − ν_L       (invariant carried from L along C⁺)
K⁺_R = θ_R + ν_R       (invariant carried from R along C⁻)

θ_P = (K⁺_R + K⁻_L) / 2
ν_P = (K⁺_R − K⁻_L) / 2
M_P = inverse_prandtl_meyer(ν_P)
μ_P = arcsin(1 / M_P)
```

**Position of P** (finite-difference approximation using average characteristic slopes):

Let:
```
s₁ = tan(θ_L + μ_L)    (slope of C⁺ from L)
s₂ = tan(θ_R − μ_R)    (slope of C⁻ from R)
```

Then:
```
x_P = (y_R − y_L + s₁·x_L − s₂·x_R) / (s₁ − s₂)
y_P = y_L + s₁·(x_P − x_L)
```

This is derived by writing the parametric equations of the two characteristic lines and solving for their intersection.

### Axis (Centerline) Point

When a C⁻ characteristic from an interior node J reaches the symmetry axis (y = 0), the boundary condition is:

```
θ = 0   (flow is parallel to the axis by symmetry)
```

The K⁺ invariant is preserved along the C⁻ characteristic from J:

```
K⁺_J = θ_J + ν_J
     = 0 + ν_axis
     → ν_axis = K⁺_J
```

Then:
```
M_axis = inverse_prandtl_meyer(ν_axis)
μ_axis = arcsin(1 / M_axis)
```

**Position of the axis node:**

```
s₂ = tan(θ_J − μ_J)           (slope of C⁻ from J)
x_axis = x_J − y_J / s₂       (extrapolate to y = 0)
y_axis = 0
```

### Wall Point

The nozzle wall is a **streamline**. By definition, the flow velocity at a streamline is tangent to the streamline — there is no flow through the wall. This means:

```
θ_W = slope angle of the wall at that point
```

The wall contour is not prescribed as an input to the solver. It is **computed** as part of the MOC solution. For a minimum-length nozzle (MLN), the wall is exactly the streamline that makes the exit flow uniform and parallel at the design Mach number.

**Flow state at a wall node W**, given the previous wall node W_prev (along the wall C⁺ characteristic) and an interior node J (on a C⁻ characteristic from the interior):

```
K⁻_prev = θ_W_prev − ν_W_prev   (from C⁺ char along wall)
K⁺_J    = θ_J + ν_J             (from C⁻ char from J)

θ_W = (K⁺_J + K⁻_prev) / 2
ν_W = (K⁺_J − K⁻_prev) / 2
M_W = inverse_prandtl_meyer(ν_W)
```

This is algebraically identical to an interior point. The difference is in the position logic:

- **First wall node** (just downstream of the throat): located at `(0, r_t)` where `r_t` is the throat radius, with `θ = θ_max` (the maximum wall angle, equal to half the total turning required for the design Mach number in a minimum-length nozzle).
- **Subsequent wall nodes**: position is the intersection of the C⁺ characteristic from W_prev and the C⁻ characteristic from J, computed with the same finite-difference formula as an interior point.

The collected sequence of wall node positions `(x_W, y_W)` defines the **nozzle wall contour** — the primary output of the MOC solver.

---

## The Characteristic Mesh Structure

For a minimum-length nozzle with n characteristic lines (n initial data points from axis to wall):

- The mesh is a **triangular grid** that contracts from the initial data line toward the exit.
- **Row 0** (initial data): n nodes from the throat wall (θ = θ_max) to the axis (θ = 0). In practice, this line is the first row of characteristics emanating from the throat.
- **Each subsequent row** has one fewer interior node than the previous, plus one axis node and one wall node.
- **Total rows**: n.
- The final row has a single node on the centerline at the nozzle exit.

Example mesh structure for n = 3:

```
Row 0 (initial):  [W0]---[3]---[2]---[1]---[A0]   (wall to axis)
                     \   / \   / \   /
Row 1:              [W1]  [P2]  [P1]  [A1]
                       \   / \   /
Row 2:                [W2]  [P3]  [A2]
                         \   /
Row 3 (exit):           [W3/Exit]---[A3]
```

Where:
- **W** = wall node (defines the nozzle contour)
- **A** = axis node
- **P** = interior node

Each row requires computing its interior nodes first, then the axis node, then the wall node. The final axis node is the exit centerline — its Mach number should equal the design exit Mach number M_e, and θ should equal 0. These two conditions serve as a check on the correctness of the computation.

The wall nodes W0 through Wn form the nozzle wall contour. Plotting these points and connecting them with a smooth curve gives the nozzle profile.

---

## For Axisymmetric Flow

The theory above applies to **2D planar flow** — an infinite slab with constant cross-section perpendicular to the x-y plane. For a real rocket nozzle, the geometry is **axisymmetric** (a body of revolution about the x-axis). In this case, the compatibility equations acquire an additional source term:

```
Along C⁺:   dθ − dν + (sin μ · sin θ) / y · ds = 0
Along C⁻:   dθ + dν − (sin μ · sin θ) / y · ds = 0
```

where `ds` is the arc-length element along the characteristic and `y` is the radial distance from the axis. This term accounts for the radial divergence of streamtubes in axisymmetric flow — as flow moves away from the axis, it occupies a larger annular area, which modifies the expansion rate.

The current project scaffolding does not include this source term. The compatibility equations in `moc/characteristics.rs` implement the 2D planar form only. **Adding axisymmetric support requires modifying the compatibility equations** to include this term and switching from a simple algebraic update to a predictor-corrector or iterative scheme along each characteristic arc.

This is documented as a **future enhancement**. The 2D planar equations are a valid starting point and produce qualitatively correct nozzle contours — the shapes will be slightly different from the true axisymmetric solution, but the flow physics and mesh structure are the same.

---

## Summary

- MOC applies because supersonic flow is governed by a **hyperbolic PDE**, which admits characteristic curves along which ODEs hold.
- There are two characteristic families: **C⁺** (left-running, carries K⁻ = θ − ν) and **C⁻** (right-running, carries K⁺ = θ + ν).
- There are three node types: **interior** (intersection of C⁺ and C⁻), **axis** (C⁻ reaches y = 0, θ = 0 by symmetry), and **wall** (streamline condition θ_W = wall slope).
- The **wall contour is not an input** — it is derived from the MOC mesh as a streamline. This is the central result: MOC tells you what wall shape to build.
- The invariant algebra in `moc/characteristics.rs` is correct. The missing pieces are:
  1. Position computation using characteristic slopes (currently the `step()` function in `moc/solver.rs` ignores slopes entirely and offsets y by a fixed 0.1 — this must be completely rewritten).
  2. Proper Mach number recovery via `gas.inverse_prandtl_meyer(nu)` after each invariant update.
- Both of those fixes depend on `prandtl_meyer()` returning radians — which requires removing the `* 180.0 / PI` conversion documented in `02_prandtl_meyer.md`.
