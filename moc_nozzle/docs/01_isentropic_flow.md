# Isentropic Flow Fundamentals

---

## What Is Isentropic Flow

**Isentropic** means two things simultaneously:

- **Adiabatic** — no heat is added to or removed from the fluid (no heat transfer through the nozzle walls).
- **Reversible** — no friction, no shock waves, no irreversible mixing. Entropy is constant.

Together these imply that the entropy of every fluid parcel is constant as it travels through the nozzle.

In reality, rocket nozzles have boundary layers, slight heat transfer, and weak oblique shocks. However, for a first-order design the isentropic assumption is excellent because:

- The core flow (away from the wall boundary layer) is nearly inviscid.
- The nozzle is designed to be shock-free (that is precisely the goal of MOC design).
- Combustion gases behave approximately as a calorically perfect gas over the temperature range of interest.

The isentropic model therefore gives accurate nozzle contours, thrust coefficients, and exit conditions for conceptual and preliminary design.

The three conservation laws still apply across any cross-section of the nozzle:

- **Continuity** — mass flow rate is constant: ρ A V = constant.
- **Momentum** — relates pressure gradient to velocity change.
- **Energy** — total enthalpy is constant (consequence of adiabatic + steady flow): h + V²/2 = h₀ = constant.

---

## Speed of Sound and Mach Number

The **speed of sound** in a calorically perfect gas is:

```
a = sqrt(γ R T)
```

where:
- `γ` is the ratio of specific heats (≈ 1.4 for air, 1.2–1.3 for hot combustion products)
- `R` is the specific gas constant (J / kg·K)
- `T` is the local static temperature (K)

Notice that `a` depends on temperature. As the gas expands and cools in the diverging section of the nozzle, the speed of sound *decreases* even as the flow velocity *increases* — which is why Mach number rises quickly.

The **Mach number** is the ratio of the local flow speed to the local speed of sound:

```
M = V / a
```

There are three regimes:

| Regime | M | Location in nozzle |
|--------|---|---------------------|
| Subsonic | M < 1 | Converging section |
| Sonic | M = 1 | Throat (minimum area) |
| Supersonic | M > 1 | Diverging section |

At M = 1 the throat is said to be **choked**. This means the mass flow rate through the nozzle has reached its maximum possible value for the given stagnation conditions — increasing the downstream pressure further cannot increase the mass flow. The MOC design operates entirely in the supersonic regime downstream of this choked throat.

---

## Stagnation (Total) Conditions

**Stagnation conditions** are the hypothetical thermodynamic state a fluid parcel would reach if it were brought to rest *isentropically* — with no heat transfer and no irreversibility. They are denoted with a subscript 0.

For a calorically perfect gas the isentropic relations give:

**Stagnation temperature:**
```
T₀/T = 1 + (γ−1)/2 · M²
```

**Stagnation pressure:**
```
P₀/P = [1 + (γ−1)/2 · M²]^(γ/(γ−1))
```

**Stagnation density:**
```
ρ₀/ρ = [1 + (γ−1)/2 · M²]^(1/(γ−1))
```

These three stagnation quantities — T₀, P₀, and ρ₀ — are **constant throughout the entire isentropic flow**. They do not change from the combustion chamber to the nozzle exit. This is extremely useful: once you know the chamber stagnation conditions and the local Mach number, you can recover the local static temperature, pressure, and density anywhere in the nozzle.

Physical interpretation: at low Mach numbers (M → 0) the static and stagnation values are nearly identical. As M increases, more of the total energy is in kinetic form, so the static temperature and pressure drop below their stagnation values.

---

## Critical (Sonic) Conditions

At the throat, M = 1. Substituting into the stagnation relations gives the **critical** (sonic) conditions, denoted with a superscript `*`:

**Critical temperature:**
```
T* = T₀ · 2/(γ+1)
```

**Critical pressure:**
```
P* = P₀ · [2/(γ+1)]^(γ/(γ−1))
```

**Critical density:**
```
ρ* = ρ₀ · [2/(γ+1)]^(1/(γ−1))
```

For air (γ = 1.4):
- T*/T₀ ≈ 0.833
- P*/P₀ ≈ 0.528
- ρ*/ρ₀ ≈ 0.634

The **throat area** A* is the reference area for the entire nozzle. Every area ratio in nozzle design is expressed relative to A*. The critical conditions are fixed for given stagnation conditions and γ — they do not depend on the nozzle geometry downstream.

---

## The Area-Mach Relation

Combining the continuity equation (ρ A V = constant) with the isentropic relations yields the **Area-Mach relation**, which connects the local cross-sectional area to the local Mach number:

```
A/A* = (1/M) · [(2/(γ+1)) · (1 + (γ−1)/2 · M²)]^((γ+1)/(2(γ−1)))
```

Key properties of this function:

- It is **monotonically decreasing** for M < 1 (the area must shrink to accelerate subsonic flow).
- It is **monotonically increasing** for M > 1 (the area must grow to accelerate supersonic flow).
- It has a **minimum of 1.0 at M = 1** (the throat).
- For any A/A* > 1 there are **two solutions**: one subsonic (M < 1) and one supersonic (M > 1).

In a converging-diverging nozzle, the converging section operates on the subsonic branch and the diverging section operates on the supersonic branch. The MOC design only needs the supersonic solution.

### Reference values for γ = 1.4

| M | A/A* |
|---|------|
| 1.5 | ≈ 1.176 |
| 2.0 | ≈ 1.687 |
| 2.5 | ≈ 2.637 |
| 3.0 | ≈ 4.235 |
| 5.0 | ≈ 25.0 |

Notice how rapidly A/A* grows with Mach number. A nozzle designed for M = 5 has an exit area 25 times the throat area.

### Inverting A/A* → M

The Area-Mach equation cannot be inverted analytically. To find M from a given A/A* on the supersonic branch, use a numerical root-finder — specifically bisection on the interval M ∈ (1, ∞). Because A/A* is strictly monotone on the supersonic branch, the root is unique and bisection converges reliably.

---

## Implementing the Area-Mach Relation in Rust

The following two methods should be added to the `GasModel` trait (and implemented for `Air`) in `core/gas.rs`.

`area_mach_ratio()` evaluates the formula directly for a given M:

```rust
/// Returns A/A* for a given Mach number (isentropic, perfect gas)
fn area_mach_ratio(&self, m: f64) -> f64 {
    let g = self.gamma();
    let t = (2.0 + (g - 1.0) * m * m) / (g + 1.0);
    (1.0 / m) * t.powf((g + 1.0) / (2.0 * (g - 1.0)))
}

/// Inverts A/A* → M (supersonic branch, M > 1)
/// Uses bisection from utils/root.rs
fn mach_from_area_ratio(&self, ae_at: f64) -> f64 {
    use crate::utils::root::bisection;
    // Area ratio is always ≥ 1.0; supersonic branch M ∈ (1, ∞)
    bisection(|m| self.area_mach_ratio(m) - ae_at, 1.0 + 1e-9, 50.0)
}
```

How `bisection` works here:

- The function passed to `bisection` is `f(M) = A/A*(M) − ae_at`.
- At the lower bracket `M = 1 + ε`: A/A* ≈ 1.0, so `f ≈ 1 − ae_at < 0` (since ae_at > 1).
- At the upper bracket `M = 50`: A/A* is enormous, so `f > 0`.
- Because A/A* is strictly increasing on the supersonic branch, the function crosses zero exactly once, and bisection will find that crossing in at most 50 halvings.
- After 50 iterations the bracket width is `(50 − 1) / 2^50 ≈ 4 × 10⁻¹⁴`, which is well within double-precision accuracy.

The upper bound of 50 is a practical ceiling — no rocket nozzle operates at M = 50 — so it is safe to use as a bracket limit.

---

## What ae_at Means for Nozzle Design

In `main.rs` the user configures the solver with `ae_at = 10.0`. This single number drives the entire design:

- The exit area is **10 times** the throat area: A_e = 10 · A*.
- Plugging ae_at = 10.0 into `mach_from_area_ratio()` with γ = 1.4 gives a design exit Mach number of approximately **M_e ≈ 3.96**.
- The entire MOC characteristic mesh is built so that the flow at the exit plane is uniform, parallel, and at exactly M_e.

This exit Mach number is the **design Mach number**. Every initial characteristic in the throat expansion fan, every interior point, and every wall boundary condition is computed to achieve this target.

Currently, `NozzleSolver` stores `ae_at` in its configuration but never calls any function to derive `M_e` from it. The `run()` method proceeds directly to the dummy solver without this critical step. Computing `M_e` via `mach_from_area_ratio()` is therefore **one of the first things to fix** when implementing the real solver.

---

## Summary

- **Isentropic flow** is adiabatic and reversible; entropy is constant. It is the standard first-order model for rocket nozzles.
- The **speed of sound** depends on temperature: a = sqrt(γ R T). Mach number M = V/a.
- **Stagnation conditions** T₀, P₀, ρ₀ are conserved throughout isentropic flow.
- The **throat** is always at M = 1 (sonic condition, choked flow). All areas are referenced to the throat area A*.
- The **Area-Mach relation** A/A*(M) uniquely determines M on each branch (subsonic or supersonic). For any A/A* > 1 there are two solutions.
- `mach_from_area_ratio()` is the first new function to implement: it inverts the Area-Mach relation on the supersonic branch using bisection, turning the configured `ae_at` into the design exit Mach number M_e that drives the rest of the MOC computation.
