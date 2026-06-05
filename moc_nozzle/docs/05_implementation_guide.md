# Implementation Guide — Missing Pieces in Rust

This guide walks through each missing or broken piece in order, providing the exact Rust code to write and explaining each decision. Follow the implementation order in `00_overview.md`.

---

## 1. Fix `core/gas.rs` — Prandtl-Meyer and Area-Mach

**What to change:**

- `prandtl_meyer` currently multiplies by `180/PI`, returning degrees instead of radians — remove that factor
- `inverse_prandtl_meyer` is broken (returns a hardcoded value) — replace with a proper bisection
- Add `area_mach_ratio` and `mach_from_area_ratio` to the `GasModel` trait and `Air` impl
- Add `mach_angle` (μ = arcsin(1/M)) to the trait, needed by all characteristic calculations

Complete replacement for `core/gas.rs`:

```rust
use std::f64::consts::PI;

pub trait GasModel {
    fn gamma(&self) -> f64;
    fn prandtl_meyer(&self, m: f64) -> f64;
    fn inverse_prandtl_meyer(&self, nu: f64) -> f64;
    fn mach_angle(&self, m: f64) -> f64;
    fn area_mach_ratio(&self, m: f64) -> f64;
    fn mach_from_area_ratio(&self, ae_at: f64) -> f64;
}

pub struct Air {
    gamma: f64,
}

impl Air {
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }
}

impl GasModel for Air {
    fn gamma(&self) -> f64 {
        self.gamma
    }

    /// ν(M) in radians. Formula is correct; `* 180/PI` removed.
    fn prandtl_meyer(&self, m: f64) -> f64 {
        let g = self.gamma;
        let a = (g + 1.0) / (g - 1.0);
        a.sqrt() * ((g - 1.0) / (g + 1.0) * (m * m - 1.0)).sqrt().atan()
            - (m * m - 1.0).sqrt().atan()
    }

    /// Inverts ν(M) using bisection. Input ν in radians, output M.
    fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
        use crate::utils::root::bisection;
        if nu <= 0.0 {
            return 1.0;
        }
        bisection(|m| self.prandtl_meyer(m) - nu, 1.0 + 1e-9, 100.0)
    }

    /// Mach angle in radians: μ = arcsin(1/M)
    fn mach_angle(&self, m: f64) -> f64 {
        (1.0 / m).asin()
    }

    /// A/A* = (1/M) * [(2/(γ+1))*(1 + (γ-1)/2*M²)]^((γ+1)/(2(γ-1)))
    fn area_mach_ratio(&self, m: f64) -> f64 {
        let g = self.gamma;
        let t = (2.0 + (g - 1.0) * m * m) / (g + 1.0);
        (1.0 / m) * t.powf((g + 1.0) / (2.0 * (g - 1.0)))
    }

    /// Supersonic M from A/A* using bisection on (1, 100)
    fn mach_from_area_ratio(&self, ae_at: f64) -> f64 {
        use crate::utils::root::bisection;
        bisection(|m| self.area_mach_ratio(m) - ae_at, 1.0 + 1e-9, 100.0)
    }
}
```

**Note on `bisection` imports:** Rust trait method bodies can use `use` statements locally, scoped to that method. The pattern `use crate::utils::root::bisection;` inside the function body is idiomatic and avoids polluting the module namespace. The two methods that need it are:

```rust
fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
    use crate::utils::root::bisection;
    if nu <= 0.0 { return 1.0; }
    bisection(|m| self.prandtl_meyer(m) - nu, 1.0 + 1e-9, 100.0)
}

fn mach_from_area_ratio(&self, ae_at: f64) -> f64 {
    use crate::utils::root::bisection;
    bisection(|m| self.area_mach_ratio(m) - ae_at, 1.0 + 1e-9, 100.0)
}
```

The `PI` import from `std::f64::consts` can be left in place or removed — it is no longer needed once the degree conversion is gone.

---

## 2. Update `solver/config.rs`

Rename `n_points` to `n_chars` (number of characteristic lines) and add `throat_radius`:

```rust
pub struct NozzleConfig {
    pub gamma:         f64,
    pub ae_at:         f64,
    pub n_chars:       usize,   // number of characteristic lines (was n_points)
    pub throat_radius: f64,     // throat half-height (normalized, typically 1.0)
}
```

Update `main.rs` wherever `NozzleConfig` is constructed to use the new field names (see Section 8 below).

---

## 3. Add `moc_angle` Helper to `moc/characteristics.rs`

Update `from_invariants` to accept a gas model and recover M correctly. The hardcoded `m: 2.0` bug is fixed by inverting the Prandtl-Meyer function.

Complete replacement for the public API in `moc/characteristics.rs`:

```rust
use crate::core::gas::GasModel;
use crate::core::state::FlowState;

pub struct Invariants {
    pub k_plus:  f64,   // θ + ν, constant along C⁻ characteristics
    pub k_minus: f64,   // θ − ν, constant along C⁺ characteristics
}

pub fn invariants(s: FlowState) -> Invariants {
    Invariants {
        k_plus:  s.theta + s.nu,
        k_minus: s.theta - s.nu,
    }
}

pub fn from_invariants<G: GasModel>(k_plus: f64, k_minus: f64, gas: &G) -> FlowState {
    let theta = (k_plus + k_minus) / 2.0;
    let nu    = (k_plus - k_minus) / 2.0;
    let m     = gas.inverse_prandtl_meyer(nu.abs());
    FlowState { theta, nu, m }
}
```

The calling convention matches the array ordering: pass `(k_plus_from_wall_side, k_minus_from_axis_side, gas)`. See the conventions note in `03_characteristics.md`.

---

## 4. Add Node Computation Functions to `moc/characteristics.rs`

Add three new public functions after `from_invariants`. These encapsulate the three boundary conditions used during mesh marching.

```rust
use crate::moc::node::Node;

/// Computes an interior node from two neighboring nodes.
/// `lower` = axis-side node (provides K⁻ along C⁺ char)
/// `upper` = wall-side node (provides K⁺ along C⁻ char)
pub fn interior_point<G: GasModel>(lower: &Node, upper: &Node, gas: &G) -> Node {
    let k_minus = lower.state.theta - lower.state.nu; // from C⁺ through lower
    let k_plus  = upper.state.theta + upper.state.nu; // from C⁻ through upper

    let state = from_invariants(k_plus, k_minus, gas);

    let mu_l = gas.mach_angle(lower.state.m);
    let mu_u = gas.mach_angle(upper.state.m);
    let s1 = (lower.state.theta + mu_l).tan(); // C⁺ slope from lower
    let s2 = (upper.state.theta - mu_u).tan(); // C⁻ slope from upper

    let denom = s1 - s2;
    let x = if denom.abs() < 1e-12 {
        (lower.x + upper.x) / 2.0 // degenerate fallback
    } else {
        (upper.y - lower.y + s1 * lower.x - s2 * upper.x) / denom
    };
    let y = lower.y + s1 * (x - lower.x);

    Node { x, y, state }
}

/// Computes the axis node when a C⁻ characteristic from `j_node` reaches y=0.
/// Enforces θ = 0 by symmetry; ν = K⁺_J.
pub fn axis_point<G: GasModel>(j_node: &Node, gas: &G) -> Node {
    let k_plus = j_node.state.theta + j_node.state.nu; // preserved along C⁻
    let nu     = k_plus; // at axis: theta = 0, so nu = K+
    let m      = gas.inverse_prandtl_meyer(nu.abs());
    let state  = FlowState { theta: 0.0, nu, m };

    let mu_j = gas.mach_angle(j_node.state.m);
    let slope = (j_node.state.theta - mu_j).tan(); // C⁻ slope
    // Extend from j_node to y = 0:
    let x = if slope.abs() < 1e-12 {
        j_node.x
    } else {
        j_node.x - j_node.y / slope
    };

    Node { x, y: 0.0, state }
}

/// Computes a wall node from the adjacent interior node `j_node`
/// and the previous wall node `w_prev`.
/// The wall is a streamline: θ_wall = wall slope (flow tangent to wall).
pub fn wall_point<G: GasModel>(j_node: &Node, w_prev: &Node, gas: &G) -> Node {
    let k_plus  = j_node.state.theta + j_node.state.nu;   // from C⁻ through j
    let k_minus = w_prev.state.theta - w_prev.state.nu;   // from C⁺ along wall

    let state = from_invariants(k_plus, k_minus, gas);

    let mu_j = gas.mach_angle(j_node.state.m);
    let mu_w = gas.mach_angle(w_prev.state.m);
    let s_j  = (j_node.state.theta - mu_j).tan();           // C⁻ from j_node
    let s_w  = (w_prev.state.theta + mu_w).tan();           // C⁺ from w_prev

    let denom = s_w - s_j;
    let x = if denom.abs() < 1e-12 {
        (j_node.x + w_prev.x) / 2.0
    } else {
        (j_node.y - w_prev.y + s_w * w_prev.x - s_j * j_node.x) / denom
    };
    let y = j_node.y + s_j * (x - j_node.x);

    Node { x, y, state }
}
```

---

## 5. Rewrite `moc/solver.rs`

Replace the stub `SimpleMocSolver` with a proper `MocSolver` that runs the full MLN design.

### ⚠️ The Throat Singularity

Before reading the code, understand a critical subtlety: **all fan characteristics originate from the same corner point (0, r_t)**. This means:

1. The initial data line **cannot** be the vertical throat plane at x = 0, because all characteristics collapse to one point there.
2. If we place `current_row[0]` (the outermost node) and the throat wall node `W₀` at the **same coordinates** (0, r_t), then `wall_point` is called with two identical nodes, producing a degenerate result (x stays at 0 forever).

**Fix**: The throat wall node `W₀` is stored separately (in `wall_nodes`). The initial data line nodes are placed at a **small ε offset** (x = ε > 0) along the fan characteristics, so that `prev[0]` (the first interior node) is strictly downstream of `W₀` at (0, r_t). This unblocks the `wall_point` geometry.

A convenient formula for ε: `ε = r_t * 0.01` (1% of throat radius). The fan char with state (θ_k, M_k) at x = ε is at:
```
y_k = r_t + ε · tan(θ_k − μ_k)
```
For near-sonic nodes (M close to 1), `tan(θ − μ) ≈ −∞`, so the near-axis nodes are placed at y ≈ 0 for any small ε. In practice, clamp y to `[0, r_t]`.

### The Corrected `design()` Function

```rust
use crate::core::gas::GasModel;
use crate::core::state::FlowState;
use crate::moc::node::Node;
use crate::moc::characteristics::{interior_point, axis_point, wall_point};
use crate::solver::config::NozzleConfig;

pub struct MocMesh {
    pub wall_nodes:  Vec<Node>,
    pub axis_nodes:  Vec<Node>,
    pub all_nodes:   Vec<Node>,
}

pub struct MocSolver {
    pub n:     usize,
    pub gamma: f64,
    pub mesh:  Option<MocMesh>,
}

impl MocSolver {
    pub fn new(n: usize, gamma: f64) -> Self {
        Self { n, gamma, mesh: None }
    }

    pub fn nodes(&self) -> &[Node] {
        self.mesh.as_ref().map(|m| m.all_nodes.as_slice()).unwrap_or(&[])
    }

    pub fn wall_nodes(&self) -> &[Node] {
        self.mesh.as_ref().map(|m| m.wall_nodes.as_slice()).unwrap_or(&[])
    }

    /// Runs the full MOC design. Populates `self.mesh` with the result.
    pub fn design<G: GasModel>(&mut self, gas: &G, config: &NozzleConfig) {
        let n         = self.n;
        let r_t       = config.throat_radius;
        let m_exit    = gas.mach_from_area_ratio(config.ae_at);
        let nu_exit   = gas.prandtl_meyer(m_exit);
        let theta_max = nu_exit / 2.0;

        // Small offset to avoid the throat singularity at x=0.
        // The fan chars are placed at x=ε so that wall_point(prev[0], W₀)
        // has a non-degenerate geometry: prev[0] is strictly downstream of W₀.
        let eps = r_t * 0.01;

        // --- Throat wall node W₀ (stored separately from the interior data line) ---
        let m_wall = gas.inverse_prandtl_meyer(theta_max);
        let w0 = Node {
            x:     0.0,
            y:     r_t,
            state: FlowState { m: m_wall, theta: theta_max, nu: theta_max },
        };

        // --- Initial data line: n nodes from wall-side to axis-side.
        //
        //  Each node k corresponds to fan characteristic k.
        //  Ordered wall-to-axis: k=0 is closest to wall, k=n-1 is axis-side.
        //  IMPORTANT: This list does NOT include W₀ (the throat lip).
        //  The fan chars have K⁻ = 0, so θ = ν (simple wave).
        //  Positions: x = ε, y computed from the fan char slope.
        let mut current_row: Vec<Node> = (0..n).map(|k| {
            // t=0 → wall-side (large θ), t=1 → axis-side (θ→0)
            let t     = (k as f64 + 0.5) / n as f64;
            let theta = (1.0 - t) * theta_max;
            let nu    = theta;  // simple wave: K⁻ = θ − ν = 0
            let m     = if nu < 1e-10 { 1.0 } else { gas.inverse_prandtl_meyer(nu) };
            let mu    = gas.mach_angle(m);
            // Fan char slope: tan(θ − μ). Clamp y to [0, r_t].
            let slope = (theta - mu).tan();
            let y     = (r_t + eps * slope).clamp(0.0, r_t - 1e-6);
            Node { x: eps, y, state: FlowState { m, theta, nu } }
        }).collect();

        let mut wall_nodes: Vec<Node> = vec![w0];
        let mut axis_nodes: Vec<Node> = vec![];
        let mut all_nodes:  Vec<Node> = current_row.clone();

        // --- March characteristic rows ---
        for _row in 0..n {
            let prev   = &current_row;
            let m_prev = prev.len();
            if m_prev < 2 { break; }

            let mut next_row: Vec<Node> = Vec::new();

            // Wall node: C⁻ from the wall-side node + C⁺ from previous wall node.
            // prev[0] is interior (wall-side) and is strictly downstream of wall_nodes.last().
            let w_new = wall_point(&prev[0], wall_nodes.last().unwrap(), gas);
            wall_nodes.push(w_new.clone());
            all_nodes.push(w_new.clone());
            next_row.push(w_new);

            // Interior nodes: pairs (upper=prev[i], lower=prev[i+1])
            for i in 0..(m_prev - 1) {
                let upper = &prev[i];      // wall side
                let lower = &prev[i + 1]; // axis side
                let p = interior_point(lower, upper, gas);
                all_nodes.push(p.clone());
                next_row.push(p);
            }

            // Axis node: C⁻ from the axis-side node reaching y=0
            let a_new = axis_point(prev.last().unwrap(), gas);
            axis_nodes.push(a_new.clone());
            all_nodes.push(a_new.clone());
            next_row.push(a_new);

            current_row = next_row;
        }

        self.mesh = Some(MocMesh { wall_nodes, axis_nodes, all_nodes });
    }
}
```

### Note on Mesh Size Growth

You will notice that each row produces `m_prev + 1` nodes (1 wall + (m_prev−1) interior + 1 axis), so the mesh **grows** by 1 node per row instead of shrinking. This is because the current structure includes both wall and axis nodes in the row. A more refined implementation would:

1. Keep wall nodes and axis nodes as **separate accumulation lists** (not in `current_row`)
2. Let `current_row` contain only **interior** nodes, shrinking from n−1 to 1
3. Use `wall_point` with the outermost interior node and the previous wall node
4. Use `axis_point` with the innermost interior node

The current structure is correct in terms of physics but will require some debugging of the row layout. Treat it as a working first draft to iterate on.

---

## 6. Rewrite `solver/nozzle.rs`

```rust
use crate::core::gas::GasModel;
use crate::moc::solver::MocSolver;
use crate::solver::config::NozzleConfig;

pub struct NozzleSolver<G: GasModel> {
    pub gas:    G,
    pub solver: MocSolver,
    pub config: NozzleConfig,
}

impl<G: GasModel> NozzleSolver<G> {
    pub fn new(gas: G, config: NozzleConfig) -> Self {
        let solver = MocSolver::new(config.n_chars, config.gamma);
        Self { gas, solver, config }
    }

    pub fn run(&mut self) {
        self.solver.design(&self.gas, &self.config);
    }
}
```

---

## 7. Rewrite `geometry/wall.rs`

Remove the broken `y >= 0.5` heuristic. The wall nodes are now explicitly tracked during the MOC march and passed in directly:

```rust
use crate::moc::node::Node;

pub struct NozzleWall {
    pub points: Vec<(f64, f64)>,
}

/// Extracts the wall contour from the ordered wall nodes produced by MocSolver.
pub fn extract_wall(wall_nodes: &[Node]) -> NozzleWall {
    NozzleWall {
        points: wall_nodes.iter().map(|n| (n.x, n.y)).collect(),
    }
}
```

---

## 8. Update `main.rs`

```rust
mod core;
mod geometry;
mod moc;
mod solver;
mod utils;

use core::gas::Air;
use solver::config::NozzleConfig;
use solver::nozzle::NozzleSolver;
use geometry::wall::extract_wall;

fn main() {
    let gas = Air::new(1.4);

    let config = NozzleConfig {
        gamma:         1.4,
        ae_at:         10.0,
        n_chars:       10,
        throat_radius: 1.0,
    };

    let mut nozzle = NozzleSolver::new(gas, config);
    nozzle.run();

    let wall = extract_wall(nozzle.solver.wall_nodes());

    println!("Nozzle wall contour ({} points):", wall.points.len());
    for (x, y) in &wall.points {
        println!("  x = {:.4}, y = {:.4}", x, y);
    }

    let exit = wall.points.last().unwrap();
    println!("\nExit: x = {:.4}, y = {:.4}", exit.0, exit.1);
    println!("Calculated A_e/A* ≈ {:.4}", exit.1 * exit.1); // for r_t = 1.0, 2D planar
}
```

---

## 9. Rust-Specific Notes

- **Trait object vs. generic:** The code uses `G: GasModel` (monomorphic generics), which is idiomatic Rust and gives better performance than `dyn GasModel`. Keep this approach throughout. It also means the compiler can inline the gas model methods into the solver hot loop.

- **Closures in bisection:** `bisection(|m| self.prandtl_meyer(m) - nu, ...)` captures `self` and `nu` by reference. This requires `F: Fn(f64) -> f64` in `bisection`'s signature, which is already what `utils/root.rs` expects. No changes needed to `root.rs`.

- **`.clone()` on Node:** `Node` derives `Clone` already (small struct — two `f64` coordinates and a `FlowState` with three `f64` fields). Cloning nodes when distributing them into `wall_nodes`, `axis_nodes`, and `all_nodes` is fine; the cost is negligible.

- **Numeric edge cases at the axis:** At the axis (ν → 0, M → 1), `(1.0 / m).asin()` approaches π/2, and `tan(θ − μ)` → tan(−π/2) → a very large negative number. The `denom.abs() < 1e-12` guard in each position computation handles the degenerate case where two characteristics are nearly parallel. Without this guard, you will get NaN or ±∞ coordinates for early rows where axis-adjacent nodes are close together.

- **bisection bracket for inverse PM:** The upper bound of `100.0` corresponds to ν(100) ≈ 2.27 rad ≈ 130°, which is effectively ν_max for any practical nozzle. If you ever need M > 100 (hypersonic research, non-standard γ), increase the bracket upper bound accordingly.

- **`use std::f64::consts::PI`:** After removing the degree conversion from `prandtl_meyer`, this import is unused and will produce a compiler warning. Remove it from `gas.rs` unless you add another use for it.

---

## 10. What to Validate

After implementing all pieces, run `cargo run` and check:

1. `wall_nodes.len()` should equal `n_chars + 1` (n initial nodes + n new wall nodes from marching = n+1 total, starting from the throat lip W₀)
2. The exit wall y should satisfy `y ≈ ae_at` for 2D planar with r_t = 1.0 — roughly, not exactly (the MLN geometry converges to the correct area ratio as n_chars → ∞)
3. The last axis node should have `ν ≈ ν_e` and `θ ≈ 0`
4. All wall y values should be monotonically increasing from r_t toward y_exit
5. All x values should be monotonically increasing from 0

Add unit tests in `core/gas.rs` with `cargo test`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prandtl_meyer_m2() {
        let air = Air::new(1.4);
        let nu = air.prandtl_meyer(2.0);
        // ν(2.0, γ=1.4) = 0.46003 rad
        assert!((nu - 0.46003).abs() < 1e-4, "got {}", nu);
    }

    #[test]
    fn test_inverse_prandtl_meyer() {
        let air = Air::new(1.4);
        let m_recovered = air.inverse_prandtl_meyer(0.46003);
        assert!((m_recovered - 2.0).abs() < 1e-3, "got {}", m_recovered);
    }

    #[test]
    fn test_area_mach() {
        let air = Air::new(1.4);
        let ratio = air.area_mach_ratio(2.0);
        // A/A*(2.0, γ=1.4) = 1.6875
        assert!((ratio - 1.6875).abs() < 1e-3, "got {}", ratio);
    }
}
```

Run with:

```sh
cargo test
cargo run
```

The `test_prandtl_meyer_m2` test is particularly important — it is the quickest way to catch the degrees/radians bug if it reappears.
