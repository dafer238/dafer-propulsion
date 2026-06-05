# The Prandtl-Meyer Function

## Physical Intuition

In supersonic flow, when a flow turns away from itself around a convex corner — an **expansion** — it accelerates. This is the opposite of what subsonic intuition suggests: in subsonic flow, acceleration requires a converging geometry (think of a narrowing pipe). In supersonic flow, the geometry opens up, and the flow fans out and speeds up to fill it.

The mechanism for this is the **Prandtl-Meyer expansion fan**: a centered expansion wave composed of infinitely many infinitesimal Mach waves. Each wave in the fan turns the flow by an infinitesimal angle dθ and simultaneously accelerates it slightly. The cumulative turning angle from M = 1 up to some Mach number M is called the **Prandtl-Meyer angle** ν(M), or the PM function.

A useful physical analogy: imagine blowing air across a curved surface that bends away from you. The air fans out as it follows the surface, and because it is spreading into a larger area, it speeds up. Each thin layer of the fan carries the flow a little further around the bend, accelerating it a little more. That is exactly what happens at a Prandtl-Meyer expansion.

This is the **fundamental mechanism inside the diverging section of a rocket nozzle**. The throat brings the flow to exactly M = 1. From that point onward, the diverging wall acts as a continuous convex surface. The nozzle wall geometry is designed so that the expansion fan terminates at the exit plane with all streamlines parallel and at the design exit Mach number M_e. The Method of Characteristics is the mathematical tool that tells you exactly what wall shape achieves this.

---

## The Prandtl-Meyer Function Formula

The Prandtl-Meyer function is defined as the total turning angle required to accelerate a flow isentropically from M = 1 to Mach number M:

```
ν(M) = √[(γ+1)/(γ−1)] · arctan(√[(γ−1)/(γ+1) · (M²−1)]) − arctan(√[M²−1])
```

### Key Properties

- **ν(1) = 0**: At sonic conditions, no turning has occurred. The flow just reached M = 1 at the throat.
- **ν increases monotonically with M**: More turning is required to reach a higher Mach number. There is a one-to-one correspondence between ν and M for M ≥ 1.
- **Maximum turning angle**: As M → ∞, ν approaches a finite maximum:

  ```
  ν_max = (π/2) · (√[(γ+1)/(γ−1)] − 1)
  ```

  For γ = 1.4: ν_max ≈ 130.5° = **2.277 rad**. No expansion fan can turn a flow by more than this, no matter how high the Mach number.

- **Units**: ν is always computed and used in **radians** for consistency with all other angle calculations throughout the codebase. This is critical — see the bug discussion below.

### Reference Table (γ = 1.4)

| M   | ν (rad) | ν (deg)         |
|-----|---------|-----------------|
| 1.0 | 0.000   | 0.00            |
| 1.5 | 0.200   | 11.50           |
| 2.0 | 0.460   | 26.38           |
| 2.5 | 0.710   | 40.68           |
| 3.0 | 0.940   | 49.76 (approx)  |
| 5.0 | 1.495   | 85.68 (approx)  |

---

## The Mach Angle

Every point in a supersonic flow has an associated **Mach angle** μ, which is the half-angle of the Mach cone — the cone within which disturbances from that point can propagate downstream:

```
μ = arcsin(1/M)
```

Key values:
- At M = 1: μ = 90° — the Mach cone is perpendicular to the flow. A normal shock is the limiting case.
- At M = 2: μ ≈ 30°
- At M = 3: μ ≈ 19.5°
- μ decreases as Mach number increases — the Mach cone narrows, and disturbances are swept more tightly downstream.

In the Method of Characteristics, the local Mach angle determines the **slope of characteristic lines** in the x-y plane. The two characteristic families have slopes `tan(θ ± μ)`, where θ is the local flow direction. This is how the geometry of the flow field is constructed from the local flow state.

---

## Current Implementation in `core/gas.rs`

### `prandtl_meyer` — Formula Correct, Units Wrong

```rust
fn prandtl_meyer(&self, m: f64) -> f64 {
    let g = self.gamma;
    let a = (g + 1.0) / (g - 1.0);
    (a.sqrt() * ( ((g - 1.0)/(g + 1.0) * (m*m - 1.0)).sqrt() ).atan()
    - (m*m - 1.0).sqrt().atan()) * 180.0 / PI
}
```

The formula itself is **mathematically correct** — it faithfully implements the Prandtl-Meyer equation. The bug is the `* 180.0 / PI` on the last line, which converts the result to **degrees**.

All other angle quantities in the codebase — `theta` (flow direction angle), characteristic slopes, and everything computed from them — are in **radians**. This mismatch breaks the Riemann invariant computations:

```
K⁺ = θ + ν
K⁻ = θ − ν
```

If θ is in radians (e.g., 0.3 rad ≈ 17°) and ν is returned in degrees (e.g., 26.38° for M = 2), then K⁺ = 0.3 + 26.38 is dimensionally incoherent and numerically wrong by roughly a factor of 57. Every downstream computation — new flow angles, new Mach numbers, characteristic slopes — will be garbage.

### `inverse_prandtl_meyer` — Broken Root-Finding

```rust
fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
    let mut m = 2.0;
    for _ in 0..50 {
        let nu_m = self.prandtl_meyer(m);
        let err = nu_m - nu;
        m -= err * 0.01;
        if m < 1.0 { m = 1.01; }
    }
    m
}
```

This function has several problems:

1. **Fixed-step gradient descent, not bisection.** The correction `err * 0.01` is an arbitrary small constant with no mathematical justification.

2. **The step size is wrong by orders of magnitude.** The derivative dν/dM ≈ 0.54 rad per unit Mach at M = 2. So for a 1-radian error, the correction is only `1.0 * 0.01 = 0.01` — moving M by just 0.01, when it needs to move by approximately `1.0 / 0.54 ≈ 1.85` units. It would take over 185 iterations just to converge by 1 unit in M, and the loop only runs 50 times.

3. **Fixed starting point at M = 2.0.** If the target ν corresponds to M = 4 or M = 5, the algorithm starts far away and takes many iterations to approach the answer — if it converges at all within 50 steps.

4. **No use of the existing bisection utility.** The `bisection` function already exists in `utils/root.rs` and is well-suited for this exactly. Since ν(M) is monotonically increasing, bisection over a bracketing interval `[1, 100]` will converge reliably in ≈ 20 iterations.

---

## Fixing the Implementation

```rust
/// Prandtl-Meyer function: returns ν(M) in **radians**
fn prandtl_meyer(&self, m: f64) -> f64 {
    let g = self.gamma;
    let a = (g + 1.0) / (g - 1.0);
    a.sqrt() * ((g - 1.0) / (g + 1.0) * (m * m - 1.0)).sqrt().atan()
        - (m * m - 1.0).sqrt().atan()
}

/// Inverse PM: returns M given ν (in radians), using bisection.
/// Valid range: ν ∈ [0, ν_max) where ν_max = (π/2)(√((γ+1)/(γ−1)) − 1)
fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
    use crate::utils::root::bisection;
    if nu <= 0.0 {
        return 1.0;
    }
    bisection(|m| self.prandtl_meyer(m) - nu, 1.0 + 1e-9, 100.0)
}
```

The changes are minimal and targeted:

- **`prandtl_meyer`**: Remove `* 180.0 / PI`. Everything else stays identical.
- **`inverse_prandtl_meyer`**: Replace the fixed-step loop with a call to the existing `bisection` utility. The bracket `[1.0 + 1e-9, 100.0]` is always valid for physical Mach numbers — `prandtl_meyer(1.0 + 1e-9)` ≈ 0 and `prandtl_meyer(100.0)` ≫ ν_max ≈ 2.277 rad, so any physical ν falls within the bracket.

The `bisection` function in `utils/root.rs` already exists and works correctly. No changes are needed there.

---

## Role in MOC: The Riemann Invariants

> The full theory is developed in [`03_moc_theory.md`](./03_moc_theory.md). This section previews only the algebraic role of ν.

In the Method of Characteristics, two quantities are conserved along the two families of characteristic lines:

- **K⁺ = θ + ν** — constant along C⁻ characteristics (right-running, "\" direction)
- **K⁻ = θ − ν** — constant along C⁺ characteristics (left-running, "/" direction)

Both θ (local flow direction angle) and ν (Prandtl-Meyer angle) must be in the **same unit** — radians — for these expressions to be dimensionally consistent and numerically correct.

Given K⁺ and K⁻ at an interior mesh point, the new flow state is recovered as:

```
θ = (K⁺ + K⁻) / 2
ν = (K⁺ − K⁻) / 2
M = inverse_prandtl_meyer(ν)
```

The Mach number is never an independent variable in MOC — it is always derived from ν via the inverse PM function. This is why fixing the units in `prandtl_meyer()` is the **first and most critical bug to fix**: every Mach number computed in the entire MOC mesh flows through these two functions.

---

## Summary

- The PM function formula in the code is mathematically correct — only the unit conversion is wrong.
- **Fix**: Remove `* 180.0 / PI` from `prandtl_meyer()` so it returns radians.
- **Fix**: Replace the fixed-step gradient descent in `inverse_prandtl_meyer()` with a call to the existing `bisection()` from `utils/root.rs`.
- After these fixes, the `from_invariants()` function in `moc/characteristics.rs` must also be updated to call `gas.inverse_prandtl_meyer(nu)` to correctly recover M from the newly computed ν at each mesh node.
