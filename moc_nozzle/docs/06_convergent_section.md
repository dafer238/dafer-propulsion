# Convergent Nozzle Section Design

## Why the Convergent Section Is Different

The Method of Characteristics (MOC) is a technique for solving **hyperbolic** partial differential equations. Supersonic flow (M > 1) is governed by a hyperbolic PDE — disturbances propagate along characteristic lines, and information travels only downstream. This is why MOC works so well for the divergent section of a nozzle.

The convergent section, however, is **subsonic** (M < 1). Subsonic flow is governed by an **elliptic** PDE — information propagates in all directions simultaneously, upstream and downstream. There are no real characteristic lines. This means:

- MOC **cannot** be applied to the convergent section
- The flow in the convergent section cannot be solved by marching downstream
- Instead, the geometry is specified directly and the flow adjusts to it

Because MOC is off the table, the convergent section is typically designed using **geometric profiles** — circular arcs, conics, or Bézier curves — chosen to be smooth and to avoid boundary layer separation. The aerodynamic shape is specified by the designer, not derived from a characteristic network.

The convergent section serves three main purposes:

1. **Transition smoothly** from the large combustion chamber (low velocity, M ≈ 0.1–0.2) to the small throat (M = 1), without abrupt changes in area that could cause separation
2. **Keep the boundary layer attached** through a smooth, gradually converging curve — sharp corners or rapid area changes would cause the boundary layer to detach
3. **Produce a uniform, one-dimensional flow at the throat** — a prerequisite for the MOC divergent section, which assumes the flow enters the divergent region as a clean, uniform sonic flow

---

## Standard Geometry: The Two-Arc Convergent Profile

The most common engineering approach uses two circular arcs connected by (optionally) a straight conical section. This gives enough geometric freedom to control the transition from the chamber to the throat while keeping the parameterization simple.

The key parameters are:

| Symbol | Meaning |
|--------|---------|
| `R_c` | Combustion chamber inner radius |
| `R_t` | Throat radius (reference; typically normalized to 1.0) |
| `R_1` | Radius of curvature of the **inlet arc** — blends the chamber wall to the cone; typically `R_1 ≈ 1.5 R_t` to `3.0 R_t` |
| `R_2` | Radius of curvature of the **throat arc** — blends the cone to the throat; typically `R_2 ≈ 0.5 R_t` to `1.5 R_t` |
| `θ_c` | Half-angle of the converging cone; typically 25° to 45° (30° is the most common choice) |

The coordinate system places `x = 0` at the throat. The divergent section occupies `x > 0`; the convergent section occupies `x < 0`.

```
 R_c ──┐
       |  arc R1
       \___
           \  straight cone θ_c
            \____
                 )──── throat at x=0, y=R_t
            arc R2
```

---

## The Three Geometric Segments

### Segment 1: Cylindrical Chamber

- Spans from `x = x_chamber` to the start of the inlet arc
- Constant radius: `y = R_c`
- Represents the combustion chamber wall — uniform cross-section

### Segment 2: Inlet Arc (tangent to chamber, tangent to cone)

- A circular arc of radius `R_1`
- **Center**: `(x_1c, R_c − R_1)` — offset inward (toward the axis) from the chamber wall
- The arc is tangent to the cylindrical chamber at one end (the tangent is horizontal, 90°) and tangent to the straight cone at the other end (at angle `90° − θ_c`)
- **Parametric form** (φ going from 0 to θ_c, measuring angle from vertical):
  - `x(φ) = x_1c + R_1 · sin(φ)`
  - `y(φ) = (R_c − R_1) + R_1 · cos(φ)`
- At φ = 0: the arc is at `(x_1c, R_c)` — tangent to the cylinder (horizontal)
- At φ = θ_c: the arc is tangent to the cone

### Segment 3: Straight Cone

- Connects the end of the inlet arc to the start of the throat arc
- Slope: `dy/dx = −tan(θ_c)` (the wall converges inward as x increases)
- This segment may have **zero length** if the inlet and throat arcs are tangent directly to each other (happens when `R_1` and `R_2` are large enough)
- The cone half-angle `θ_c` controls how steeply the nozzle converges

### Segment 4: Throat Arc (tangent to cone, tangent to throat)

- A circular arc of radius `R_2`
- **Center**: `(0, R_t + R_2)` — directly above the throat on the axis
- The arc is tangent to the straight cone at one end and ends **horizontally** at the throat
- **Parametric form** (φ going from 0 to θ_c, starting at the throat and moving upstream):
  - `x(φ) = −R_2 · sin(φ)`
  - `y(φ) = R_t + R_2 − R_2 · cos(φ)`
- At φ = 0: the arc is at `(0, R_t)` — the throat, horizontal tangent
- At φ = θ_c: the arc is tangent to the straight cone

> **Important**: the throat arc ends exactly at `(0, R_t)` with a horizontal tangent. This point is the start of the divergent section and becomes the first wall node in the MOC solution. The horizontal tangent ensures continuity of both position and slope across the throat.

---

## Parameter Selection Guidelines

| Parameter | Conservative | Typical | Aggressive |
|-----------|-------------|---------|------------|
| θ_c | 20° | 30° | 45° |
| R_1 / R_t | 3.0 | 1.5 | 0.8 |
| R_2 / R_t | 1.5 | 0.8 | 0.5 |
| R_c / R_t | 3.0 | 2.5 | 1.5 |

**Rules of thumb:**

- **Larger R_1 and R_2** → smoother curvature, lower risk of flow separation, but longer convergent section
- **Smaller θ_c** → more gradual convergence, longer nozzle overall; large θ_c (> 45°) risks separation
- **For initial design**, the following are reliable starting values:
  - `θ_c = 30°`
  - `R_1 = 1.5 · R_t`
  - `R_2 = 1.5 · R_t`
  - `R_c = 2.5 · R_t`

These values produce a compact nozzle with smooth transitions that works well across a wide range of pressure ratios and chamber conditions.

---

## Rust Implementation

Create a new file `src/geometry/convergent.rs` with the following complete implementation:

```rust
use std::f64::consts::PI;

/// Parameters for a two-arc convergent nozzle section.
pub struct ConvergentConfig {
    /// Throat radius (normalized, typically 1.0)
    pub r_throat: f64,
    /// Combustion chamber inner radius
    pub r_chamber: f64,
    /// Radius of curvature of inlet arc (blends chamber to cone)
    pub r_inlet_arc: f64,
    /// Radius of curvature of throat arc (blends cone to throat)
    pub r_throat_arc: f64,
    /// Cone half-angle in radians (e.g., 30° → PI/6)
    pub cone_half_angle: f64,
    /// Number of points per segment for discretization
    pub n_points: usize,
}

impl ConvergentConfig {
    /// Typical rocket nozzle convergent section defaults
    pub fn typical(r_throat: f64) -> Self {
        Self {
            r_throat,
            r_chamber:       2.5 * r_throat,
            r_inlet_arc:     1.5 * r_throat,
            r_throat_arc:    1.5 * r_throat,
            cone_half_angle: 30.0_f64.to_radians(),
            n_points:        20,
        }
    }
}

/// A 2D contour point (axial position x, radial position y)
#[derive(Clone, Debug)]
pub struct ContourPoint {
    pub x: f64,
    pub y: f64,
}

/// Generates the convergent nozzle wall contour (from chamber to throat).
///
/// Returns points ordered from chamber (most negative x) to throat (x = 0).
/// The throat is at (0, r_throat), matching the divergent section start.
///
/// Coordinate system: x = 0 at throat, x < 0 in the convergent section.
pub fn convergent_contour(cfg: &ConvergentConfig) -> Vec<ContourPoint> {
    let theta_c = cfg.cone_half_angle;
    let r_t = cfg.r_throat;
    let r_c = cfg.r_chamber;
    let r1  = cfg.r_inlet_arc;
    let r2  = cfg.r_throat_arc;
    let n   = cfg.n_points;

    let mut points = Vec::new();

    // ── Segment 4: Throat arc (closest to throat, x near 0) ──────────────
    // Center at (0, r_t + r2). Arc from angle (π/2 + θ_c) down to π/2.
    // At φ=0 (angle = π/2 + θ_c from center): tangent to cone
    // At φ=θ_c (angle = π/2 from center): horizontal, at (0, r_t)
    let throat_arc_center_x = 0.0_f64;
    let throat_arc_center_y = r_t + r2;

    // Sample from φ=θ_c down to φ=0 (reversed so we output chamber→throat)
    // We'll reverse the whole contour at the end anyway — build throat→chamber first
    // then reverse.

    // Build in order: throat-end first (x=0), then going upstream.
    // Throat arc: φ from 0 to θ_c means angle from (π/2) to (π/2 + θ_c)
    for i in 0..=n {
        let t   = i as f64 / n as f64;
        let phi = t * theta_c;
        let x = throat_arc_center_x - r2 * phi.sin();
        let y = throat_arc_center_y - r2 * phi.cos();
        points.push(ContourPoint { x, y });
    }

    // Point where throat arc ends (junction with cone):
    let x_throat_arc_end = -r2 * theta_c.sin();
    let y_throat_arc_end =  throat_arc_center_y - r2 * theta_c.cos();

    // ── Segment 3: Straight cone ──────────────────────────────────────────
    // Slope: dy/dx = -tan(θ_c) going upstream (x decreasing, y increasing)
    // End of cone (at inlet arc junction):
    // The inlet arc has center at (x1c, r_c - r1).
    // Arc ends tangent to cone at angle (π/2 - θ_c) from center.
    let x_inlet_arc_center_y_part = r_c - r1;
    // The tangent point on inlet arc: 
    //   x_tp = x1c + r1 * cos(π/2 - θ_c) = x1c + r1*sin(θ_c)  [negative because left of center]
    //   y_tp = (r_c - r1) + r1 * sin(π/2 - θ_c) = (r_c - r1) + r1*cos(θ_c)
    // But we need to find x1c first.
    // The straight cone passes through (x_throat_arc_end, y_throat_arc_end)
    // with slope tan(θ_c) (in the y-direction, dy/dx = tan(θ_c) going upstream, 
    // i.e., as x decreases, y increases).
    // Cone line: y - y_throat_arc_end = tan(θ_c) * (x_throat_arc_end - x)   [x < x_throat_arc_end is further upstream]
    // Equivalently: y = y_throat_arc_end + tan(θ_c) * (x_throat_arc_end - x)

    // Inlet arc tangent point y:
    let y_inlet_tangent = x_inlet_arc_center_y_part + r1 * theta_c.cos();
    // From cone equation:
    let x_inlet_tangent = x_throat_arc_end - (y_inlet_tangent - y_throat_arc_end) / theta_c.tan();

    // Inlet arc center x:
    // Center is offset perpendicular to cone at the tangent point.
    // The perpendicular to cone (slope tan θ_c) has slope -1/tan(θ_c) = -cot(θ_c)
    // Center = tangent_point + r1 * (unit normal toward axis side)
    // Normal toward axis (downward-left): direction (-sin θ_c, -cos θ_c)... 
    // Actually: the normal pointing toward lower y is (-sin θ_c, cos θ_c) rotated...
    // Easier: we know center_y = r_c - r1, and center_x = x_inlet_tangent - r1*sin(θ_c)
    // (since from center to tangent point: direction = (sin θ_c, cos θ_c) * r1 at this angle)
    let x_inlet_arc_center = x_inlet_tangent - r1 * theta_c.sin();

    // Cone segment: from (x_throat_arc_end, y_throat_arc_end) to (x_inlet_tangent, y_inlet_tangent)
    // Parameterize by y going from y_throat_arc_end to y_inlet_tangent
    let cone_len = (x_inlet_tangent - x_throat_arc_end).hypot(y_inlet_tangent - y_throat_arc_end);
    if cone_len > 1e-6 {
        for i in 1..=n {
            let t = i as f64 / n as f64;
            let y = y_throat_arc_end + t * (y_inlet_tangent - y_throat_arc_end);
            let x = x_throat_arc_end - (y - y_throat_arc_end) / theta_c.tan();
            points.push(ContourPoint { x, y });
        }
    }

    // ── Segment 2: Inlet arc ──────────────────────────────────────────────
    // Center at (x_inlet_arc_center, r_c - r1).
    // Arc from angle (π/2 - θ_c) [tangent to cone] to π/2 [tangent to cylinder].
    for i in 1..=n {
        let t   = i as f64 / n as f64;
        let phi = theta_c * (1.0 - t); // from θ_c down to 0
        let x   = x_inlet_arc_center + r1 * phi.sin();
        let y   = x_inlet_arc_center_y_part + r1 * phi.cos();
        points.push(ContourPoint { x, y });
    }

    // ── Segment 1: Cylindrical chamber (short stub upstream) ──────────────
    // Start of arc 1: at (x_inlet_arc_center, r_c) — horizontal tangent
    let x_chamber_start = x_inlet_arc_center - 1.5 * r_t; // a short chamber section
    let x_chamber_end   = x_inlet_arc_center;
    for i in 1..=n {
        let t = i as f64 / n as f64;
        let x = x_chamber_end + (x_chamber_start - x_chamber_end) * t;
        points.push(ContourPoint { x, y: r_c });
    }

    // Reverse so output is ordered chamber → throat (increasing x toward 0)
    points.reverse();
    points
}
```

---

## Adding to `geometry/mod.rs`

Update `src/geometry/mod.rs` to expose the new module:

```rust
pub mod wall;
pub mod convergent;
```

---

## Connecting to the Full Nozzle

With both sections implemented, the full nozzle profile is assembled by concatenating the two contours:

1. **Generate the convergent contour** — points ordered from `x_chamber` (far left) to `x = 0`, `y = R_t`
2. **Generate the divergent contour** — wall nodes from the MOC solution, starting at `x = 0`, `y = R_t`
3. **Concatenate** — convergent contour followed by divergent contour gives the complete wall profile

The critical joining condition: both sections must meet at `(0.0, r_t)` with θ = 0 (a horizontal tangent) at the throat. The throat arc is constructed to end horizontally by design, and the first wall node of the divergent MOC section also lies at the throat with θ ≈ 0. This guarantees a smooth, continuous wall with no kink at the junction.

Here is a snippet for assembling and printing the full profile in `main.rs`:

```rust
use geometry::convergent::{convergent_contour, ConvergentConfig};
use geometry::wall::extract_wall;

// Divergent section (from MOC)
let div_wall = extract_wall(nozzle.solver.wall_nodes());

// Convergent section
let conv_cfg = ConvergentConfig::typical(config.throat_radius);
let conv_wall = convergent_contour(&conv_cfg);

// Full nozzle: convergent (x < 0) + divergent (x ≥ 0)
println!("=== Convergent Section ===");
for p in &conv_wall {
    println!("  x = {:.4},  y = {:.4}", p.x, p.y);
}
println!("=== Divergent Section ===");
for (x, y) in &div_wall.points {
    println!("  x = {:.4},  y = {:.4}", x, y);
}
```

---

## Physical Meaning: What the Convergent Section Does

Understanding the physics helps build intuition for why the geometry choices matter:

- The **combustion chamber** contains high-pressure gas at low velocity (M ≈ 0.1–0.2). The large cross-sectional area means the gas barely moves.
- As the cross-section shrinks through the converging section, the gas accelerates. The **isentropic flow relations** connect area ratio to Mach number — for subsonic flow, smaller area means higher Mach number.
- The **throat** is the narrowest point. By the isentropic area-Mach relation, the Mach number reaches exactly **M = 1** at the throat. This is called **choked flow**.
- The **mass flow rate** through the nozzle is set entirely by the throat area: `ṁ = ρ* · a* · A*`, where the starred quantities are the sonic (throat) conditions. The throat is the bottleneck.
- Once the flow is choked, **downstream conditions cannot propagate upstream**. Changes in back pressure, altitude, or exit geometry have no effect on the upstream combustion chamber or the throat mass flow. This is why rocket engines are largely insensitive to altitude changes in thrust once the nozzle is choked — the throat fixes the flow.
- The **smooth throat arc** (governed by `R_2`) is critical: a sharp corner at the throat would generate an oblique shock immediately downstream, disrupting the clean supersonic flow that the MOC divergent section assumes. The radius `R_2` softens this corner and allows the flow to turn smoothly from subsonic convergence to supersonic expansion.

---

## Rao's Thrust-Optimized Nozzle vs. MLN

The divergent section of this project is designed as a **Minimum-Length Nozzle (MLN)**. It is worth understanding how this compares to the more commonly used bell nozzle:

**Minimum-Length Nozzle (MLN)**
- Produces the **shortest possible nozzle** for a given exit Mach number `M_e`
- Uses a **sharp corner** at the throat wall with a wall angle `θ_max = ν(M_e) / 2`, where `ν` is the Prandtl–Meyer function
- The sharp corner generates a centered Prandtl–Meyer fan that fills the nozzle
- Produces **perfectly uniform, parallel exit flow** (zero flow divergence losses)
- Drawback: the short length means a very aggressive expansion; the sharp throat corner is geometrically idealized and not physically achievable with zero radius

**Rao's Thrust-Optimized Parabolic Nozzle (Bell Nozzle)**
- **Longer** than the MLN but produces more total thrust by optimizing the trade-off between exit flow uniformity and nozzle length
- Uses a **parabolic approximation** to the wall contour derived from an optimization that maximizes the thrust integral subject to exit flow constraints
- The familiar "bell nozzle" shape — the standard for most liquid rocket engines
- **80% bell** and **60% bell** variants are common (percentage refers to length relative to an equivalent 15° conical nozzle)
- The optimization requires solving for the wall shape that maximizes `∫ p dA` over the exit plane, which is significantly more mathematically involved than MOC alone

The current project implements the MLN. Extending it to Rao's nozzle would require adding an optimization loop over the exit boundary conditions, which is a substantial additional step.

---

## Summary

- The convergent section is **subsonic** — the MOC does not apply; there are no real characteristics
- The standard engineering approach uses a **two-arc model**: an inlet arc (blending chamber to cone) and a throat arc (blending cone to throat), connected by a straight conical section
- **Key parameters**: `R_c`, `R_t`, `R_1`, `R_2`, `θ_c`
- **Typical values**: `θ_c = 30°`, `R_1 = R_2 = 1.5 R_t`, `R_c = 2.5 R_t`
- The throat arc is constructed to end at `(0, R_t)` with a **horizontal tangent** — this is the initial wall point for the divergent MOC section, ensuring a smooth junction
- For a complete nozzle profile: **concatenate convergent + divergent contours** at the throat point
- The current project focuses on the divergent section; the `geometry/convergent.rs` module described here adds the full profile capability
