# V2 — Variable Specific Heats via NASA Polynomials

## Why γ Is Not Constant

In an introductory thermodynamics course, the specific heat ratio γ = Cp/Cv is treated as a constant
(1.4 for air at room temperature). This is a **calorically perfect gas** assumption — it holds only
when the gas molecules have no vibrational or rotational energy modes excited.

For real rocket propellants:
- Combustion produces hot triatomic molecules (H₂O, CO₂) which have active vibrational modes at
  high temperature
- γ for H₂O drops from 1.33 at 1000 K to ~1.15 at 3000 K
- The gas cools from ~3000 K (chamber) to ~500–1000 K (exit) during expansion
- Using γ = const = 1.2 (a typical "frozen" average) is better than 1.4, but still wrong at every
  point

A **thermally perfect gas** has:
- Cp = Cp(T) (temperature-dependent)
- Perfect gas equation of state still holds: P = ρ R T / M_mol
- Speed of sound: a² = γ(T) · R · T / M_mol

The difference in calculated exit Mach for ae_at = 10 between γ = 1.4 (air) and γ = 1.2 (typical
rocket gas) can exceed 15%. The nozzle contour shape changes significantly.

---

## The NASA-7 Polynomial Format

NASA maintains a database of thermodynamic properties for hundreds of gas-phase species, stored as
polynomial fits to Cp(T), h(T), and s(T).

**The NASA-7 format** (7 coefficients per temperature range):

For a temperature range [T_low, T_high]:

```
Cp(T)/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴
h(T)/(RT) = a₁ + a₂T/2 + a₃T²/3 + a₄T³/4 + a₅T⁴/5 + a₆/T
s°(T)/R   = a₁·ln(T) + a₂T + a₃T²/2 + a₄T³/3 + a₅T⁴/4 + a₇
```

Where R = 8.314 J/(mol·K) is the universal gas constant.

Most species have two temperature ranges (typically 200–1000 K and 1000–6000 K). The two sets of
7 coefficients ensure continuity at the boundary temperature.

From the coefficients, derive γ(T):

```
Cv(T)/R = Cp(T)/R - 1   (for an ideal gas: Cp - Cv = R)
γ(T) = Cp(T)/Cv(T) = Cp(T) / (Cp(T) - R)
```

**Example: H₂O coefficients** (1000–6000 K, from NASA TM-4513):

```
a₁ = 2.67703787
a₂ = 2.97318329e-3
a₃ = -7.73769690e-7
a₄ = 9.44334653e-11
a₅ = -4.26900959e-15
a₆ = -2.98858938e4
a₇ = 6.88255571
```

For H₂O at T = 3000 K:
- Cp/R ≈ 2.677 + 2.973e-3 × 3000 + ... ≈ 6.55
- γ = Cp/(Cp - 1) = 6.55/5.55 ≈ 1.18

**Where to find NASA polynomial coefficients:**
- NASA CEA database: bundled with NASA CEA executable (free, downloadable)
- NIST Chemistry WebBook: https://webbook.nist.gov
- Burcat database (Burcat & Ruscic): http://garfield.chem.elte.hu/Burcat/NEWNASA.OLD
- The Python package `cantera` contains the full JANAF database

---

## How γ(T) Affects Nozzle Flow

In isentropic flow with variable γ, the stagnation relations become more complex.

The energy equation along a streamline is still:

```
h(T₀) = h(T) + V²/2     (h = enthalpy, V = flow speed)
```

But now h is not simply Cp·T — it is the integral:

```
h(T) = h(T_ref) + ∫_{T_ref}^{T} Cp(T') dT'
```

The isentropic condition P·ρ^{-γ} = const is **NOT valid** when γ = γ(T). Instead, isentropy is
expressed through the entropy differential:

```
ds = Cp(T) dT/T - R dP/P = 0   (isentropic → ds = 0)
```

The **T-P relation** for isentropic expansion with variable γ:

```
∫_{T₀}^{T} Cp(T')/T' dT' = R · ln(P/P₀)
```

This must be integrated numerically. The familiar formula T/T₀ = (P/P₀)^{(γ-1)/γ} is only the
constant-γ special case of this integral.

Similarly, the **T-M relation** comes from combining:
1. Energy: h(T₀) - h(T) = V²/2 = M²·a²(T)/2 = M²·γ(T)·R·T / (2·M_mol)
2. Definition of Mach: M = V/a(T) where a(T) = √(γ(T)·R·T/M_mol)

These two equations together determine T as a function of M for a given stagnation state — but
the relationship is implicit and must be solved numerically.

---

## The Generalized Prandtl-Meyer Function

For constant γ, ν(M) has the closed-form formula implemented in `core/gas.rs`. For variable γ,
the formula becomes an integral over Mach number:

```
ν(M₁ → M₂) = ∫_{M₁}^{M₂} √(M²-1) / (M · (1 + (γ(T(M))-1)/2 · M²)) dM
```

where T(M) is found from the energy equation (requires solving for T at each M value in the
integrand). The full evaluation chain for a single PM function call is therefore:

1. For each quadrature point Mᵢ, call the T(M) solver → Tᵢ
2. Evaluate γ(Tᵢ) via NASA polynomials
3. Evaluate the integrand f(Mᵢ)
4. Accumulate the weighted sum (e.g., Simpson's rule)

**Numerical evaluation** (Simpson's rule with N = 100 points):

```
ν(1 → M) ≈ Σᵢ wᵢ · f(Mᵢ)
where f(M) = √(M²-1) / (M · (1 + (γ(T(M))-1)/2 · M²))
```

> **Note:** For γ slowly varying (changes less than ~10% across the expansion), an approximation
> using a mean effective γ_eff is acceptable:
>
> ```
> γ_eff = average of γ(T(M)) over the expected M range
> ```
>
> This simplifies implementation significantly: use the existing PM formula with γ = γ_eff instead
> of γ = 1.4. This is the recommended starting point for V2 before committing to the full numerical
> integration.

---

## Computing T(M) from the Energy Equation

Given stagnation temperature T₀ and a local Mach number M, we want to find the static temperature T.

The energy equation requires:

```
h(T₀) - h(T) = M²·γ(T)·R·T / (2·M_mol)
```

This is a **nonlinear equation in T** because:
- γ(T) = Cp(T)/(Cp(T) - R/M_mol) depends on T through the NASA polynomial
- h(T) = ∫Cp dT also depends on T through the NASA polynomial

We solve it with bisection. For a known T₀ and M, the algorithm is:

1. Compute h(T₀) using NASA polynomials
2. Define the residual function:
   ```
   g(T) = h(T₀) - h(T) - M²·γ(T)·R·T / (2·M_mol)
   ```
3. Bracket: T ∈ [T_min, T₀]
   - At T = T₀: g(T₀) = 0 - M²·(...) < 0 (for M > 0)
   - At T = T_min: h(T₀) - h(T_min) is large and positive, so g > 0
4. Apply bisection to find the root

The bracket always exists because the left-hand side h(T₀) - h(T) is monotonically increasing
as T decreases from T₀, while the right-hand side M²·γ(T)·R·T/(2·M_mol) goes to zero as T → 0.
The crossing is guaranteed.

**Practical temperature bounds:**
- T_min ≈ 200 K (below this, NASA polynomial fits become unreliable and the gas condenses)
- T_max = T₀ (upper bound, corresponding to M = 0)

---

## Rust Implementation

### New trait: `ThermodynamicGas`

This trait extends the base `GasModel` with temperature-dependent property evaluation. Placed in
`src/core/thermo.rs`:

```rust
/// Extension of GasModel for temperature-dependent properties.
///
/// Implementors provide Cp(T) via NASA polynomials; this trait then derives
/// γ(T), sound speed, and enthalpy from those building blocks.
pub trait ThermodynamicGas: GasModel {
    /// Cp(T)/R — dimensionless specific heat at constant pressure.
    fn cp_over_r(&self, t: f64) -> f64;

    /// Enthalpy h(T) in J/kg.
    fn enthalpy(&self, t: f64) -> f64;

    /// Local γ(T) = Cp(T) / (Cp(T) - R/M_mol).
    ///
    /// Note: cp_over_r returns Cp divided by R_universal (8.314 J/mol·K),
    /// so the ideal-gas relation Cp - Cv = R gives Cv/R = Cp/R - 1 per mole.
    fn gamma_at_t(&self, t: f64) -> f64 {
        let cp_r = self.cp_over_r(t);
        cp_r / (cp_r - 1.0)
    }

    /// Speed of sound at temperature T in m/s.
    fn sound_speed(&self, t: f64) -> f64 {
        let r_specific = 8314.0 / self.molar_mass();  // J/(kg·K)
        (self.gamma_at_t(t) * r_specific * t).sqrt()
    }

    /// Molar mass of the gas (or mixture) in g/mol.
    fn molar_mass(&self) -> f64;

    /// Compute the stagnation temperature T₀ given static temperature T and Mach M.
    ///
    /// Solves: h(T₀) = h(T) + M²·γ(T)·R·T / 2  for T₀.
    fn stagnation_temperature(&self, t: f64, m: f64) -> f64 {
        let r_sp = 8314.0 / self.molar_mass();
        let ke = 0.5 * m * m * self.gamma_at_t(t) * r_sp * t;
        // h(T₀) = h(T) + ke → bisect to find T₀
        use crate::utils::root::bisection;
        let h_target = self.enthalpy(t) + ke;
        bisection(|t0| self.enthalpy(t0) - h_target, t, t * 5.0)
    }

    /// Compute the static temperature T given stagnation temperature T₀ and Mach M.
    ///
    /// Solves: h(T₀) - h(T) = M²·γ(T)·R·T / 2  for T using bisection.
    fn static_temperature(&self, t0: f64, m: f64, t_min: f64) -> f64 {
        let r_sp = 8314.0 / self.molar_mass();
        let h0 = self.enthalpy(t0);
        use crate::utils::root::bisection;
        bisection(
            |t| h0 - self.enthalpy(t) - 0.5 * m * m * self.gamma_at_t(t) * r_sp * t,
            t_min,
            t0,
        )
    }

    /// Numerically integrate the Prandtl-Meyer function from M=1 to M using
    /// the temperature-dependent integrand.
    ///
    /// Uses Simpson's rule with N steps. Requires `static_temperature` for
    /// each quadrature point, so this is O(N * bisection_iters) per call.
    fn prandtl_meyer_numerical(&self, m: f64, t0: f64, t_min: f64, n: usize) -> f64 {
        assert!(m >= 1.0, "PM function undefined for M < 1");
        if (m - 1.0).abs() < 1e-12 { return 0.0; }

        let m_start = 1.0 + 1e-9;  // avoid singularity at M=1
        let dm = (m - m_start) / n as f64;

        let integrand = |mi: f64| -> f64 {
            let ti = self.static_temperature(t0, mi, t_min);
            let gi = self.gamma_at_t(ti);
            let num = (mi * mi - 1.0).sqrt();
            let den = mi * (1.0 + 0.5 * (gi - 1.0) * mi * mi);
            num / den
        };

        // Simpson's rule (n must be even; silently use n+1 if odd)
        let steps = if n % 2 == 0 { n } else { n + 1 };
        let mut sum = integrand(m_start) + integrand(m_start + steps as f64 * dm);
        for i in 1..steps {
            let mi = m_start + i as f64 * dm;
            let weight = if i % 2 == 0 { 2.0 } else { 4.0 };
            sum += weight * integrand(mi);
        }
        sum * dm / 3.0
    }
}
```

### New struct: `NasaPolynomial`

The lowest-level building block. Placed in `src/core/nasa.rs`:

```rust
/// Represents a single NASA-7 polynomial valid for one temperature range.
///
/// The 7 coefficients encode Cp(T)/R, h(T)/(RT), and s°(T)/R as polynomial
/// fits. Most species have two of these (low-T and high-T ranges) stored
/// together in a `NasaSpecies`.
#[derive(Clone, Debug)]
pub struct NasaPolynomial {
    pub t_low:  f64,       // lower bound of temperature range [K]
    pub t_high: f64,       // upper bound of temperature range [K]
    pub a:      [f64; 7],  // coefficients a1..a7
}

impl NasaPolynomial {
    /// Cp(T)/R using the polynomial fit.
    ///
    /// Dimensionless. Multiply by R_universal (8.314 J/mol·K) for Cp in J/(mol·K),
    /// or by R_specific = R_universal / M_mol for Cp in J/(kg·K).
    pub fn cp_over_r(&self, t: f64) -> f64 {
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;
        self.a[0]
            + self.a[1] * t
            + self.a[2] * t2
            + self.a[3] * t3
            + self.a[4] * t4
    }

    /// h(T)/(R·T) — dimensionless reduced enthalpy.
    ///
    /// Includes the integration constant a₆/T which encodes the heat of
    /// formation at 298 K. To get h in J/mol: multiply by R * T.
    pub fn h_over_rt(&self, t: f64) -> f64 {
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;
        self.a[0]
            + self.a[1] * t / 2.0
            + self.a[2] * t2 / 3.0
            + self.a[3] * t3 / 4.0
            + self.a[4] * t4 / 5.0
            + self.a[5] / t
    }

    /// Standard entropy s°(T)/R — dimensionless.
    ///
    /// Reference state is 1 atm. Useful for Gibbs free energy calculations
    /// in equilibrium chemistry (V4+). Not needed for frozen-flow MOC.
    pub fn s_over_r(&self, t: f64) -> f64 {
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;
        self.a[0] * t.ln()
            + self.a[1] * t
            + self.a[2] * t2 / 2.0
            + self.a[3] * t3 / 3.0
            + self.a[4] * t4 / 4.0
            + self.a[6]
    }
}
```

### New struct: `NasaSpecies`

Wraps one or two `NasaPolynomial` instances for a single chemical species, with automatic range
selection. Also placed in `src/core/nasa.rs`:

```rust
/// A single chemical species with NASA-7 polynomial thermodynamic data.
///
/// Stores one polynomial per temperature range (typically two: 200–1000 K
/// and 1000–6000 K). The ranges should be contiguous and non-overlapping.
#[derive(Clone, Debug)]
pub struct NasaSpecies {
    pub name:        String,
    pub molar_mass:  f64,                    // g/mol
    pub polynomials: Vec<NasaPolynomial>,    // sorted ascending by t_low
}

impl NasaSpecies {
    /// Select the correct polynomial for temperature T.
    ///
    /// Clamps to the nearest range if T falls outside all defined ranges,
    /// rather than panicking — necessary because edge temperatures near
    /// 200 K or 6000 K may be just outside bounds due to floating-point.
    fn poly_at(&self, t: f64) -> &NasaPolynomial {
        for p in &self.polynomials {
            if t >= p.t_low && t <= p.t_high {
                return p;
            }
        }
        // Clamp to nearest range
        if t < self.polynomials[0].t_low {
            &self.polynomials[0]
        } else {
            self.polynomials.last().unwrap()
        }
    }

    /// Dimensionless Cp(T)/R for this species.
    pub fn cp_over_r(&self, t: f64) -> f64 {
        self.poly_at(t).cp_over_r(t)
    }

    /// Dimensionless h(T)/(R·T) for this species.
    pub fn h_over_rt(&self, t: f64) -> f64 {
        self.poly_at(t).h_over_rt(t)
    }

    /// Absolute enthalpy h(T) in J/mol.
    ///
    /// Includes heat of formation. For relative enthalpy differences
    /// (which is all the MOC energy equation needs), the absolute
    /// reference cancels out as long as you use the same species data
    /// at both T₀ and T.
    pub fn enthalpy_j_per_mol(&self, t: f64) -> f64 {
        self.h_over_rt(t) * 8.314 * t
    }

    /// Absolute enthalpy h(T) in J/kg.
    pub fn enthalpy_j_per_kg(&self, t: f64) -> f64 {
        self.enthalpy_j_per_mol(t) / (self.molar_mass * 1e-3)
    }

    /// Local γ(T) for this species alone.
    pub fn gamma_at_t(&self, t: f64) -> f64 {
        let cp_r = self.cp_over_r(t);
        cp_r / (cp_r - 1.0)
    }
}
```

### Hard-coded reference data for common rocket species

For V2, include coefficients for the five most common LOX/LH₂ and LOX/RP-1 combustion products
directly in `src/core/nasa.rs`. These save users from needing to locate and parse the CEA database
files during early development:

```rust
/// Returns NASA-7 polynomial data for H₂O (water vapor).
/// Source: NASA TM-4513 (McBride, Gordon, Reno 1993).
pub fn h2o() -> NasaSpecies {
    NasaSpecies {
        name: "H2O".into(),
        molar_mass: 18.015,
        polynomials: vec![
            // Low-T range: 200–1000 K
            NasaPolynomial {
                t_low: 200.0, t_high: 1000.0,
                a: [4.19864056e0, -2.03643410e-3, 6.52040211e-6,
                    -5.48797062e-9, 1.77197250e-12, -3.02937267e4, -8.49032208e-1],
            },
            // High-T range: 1000–6000 K
            NasaPolynomial {
                t_low: 1000.0, t_high: 6000.0,
                a: [2.67703787e0, 2.97318329e-3, -7.73769690e-7,
                    9.44334653e-11, -4.26900959e-15, -2.98858938e4, 6.88255571e0],
            },
        ],
    }
}

/// Returns NASA-7 polynomial data for H₂ (hydrogen gas).
pub fn h2() -> NasaSpecies {
    NasaSpecies {
        name: "H2".into(),
        molar_mass: 2.016,
        polynomials: vec![
            NasaPolynomial {
                t_low: 200.0, t_high: 1000.0,
                a: [2.34433112e0, 7.98052075e-3, -1.94781510e-5,
                    2.01572094e-8, -7.37611761e-12, -9.17935173e2, 6.83010238e-1],
            },
            NasaPolynomial {
                t_low: 1000.0, t_high: 6000.0,
                a: [2.93286575e0, 8.26607967e-4, -1.46402364e-7,
                    1.54100414e-11, -6.88804800e-16, -8.13065597e2, -1.02432865e0],
            },
        ],
    }
}

/// Returns NASA-7 polynomial data for CO₂ (carbon dioxide).
pub fn co2() -> NasaSpecies {
    NasaSpecies {
        name: "CO2".into(),
        molar_mass: 44.010,
        polynomials: vec![
            NasaPolynomial {
                t_low: 200.0, t_high: 1000.0,
                a: [2.35677352e0, 8.98459677e-3, -7.12356269e-6,
                    2.45919022e-9, -1.43699548e-13, -4.83719697e4, 9.90105222e0],
            },
            NasaPolynomial {
                t_low: 1000.0, t_high: 6000.0,
                a: [4.63659493e0, 2.74146090e-3, -9.95897590e-7,
                    1.60391440e-10, -9.16198570e-15, -4.90249341e4, -1.93534855e0],
            },
        ],
    }
}

/// Returns NASA-7 polynomial data for CO (carbon monoxide).
pub fn co() -> NasaSpecies {
    NasaSpecies {
        name: "CO".into(),
        molar_mass: 28.010,
        polynomials: vec![
            NasaPolynomial {
                t_low: 200.0, t_high: 1000.0,
                a: [3.57953347e0, -6.10353680e-4, 1.01681433e-6,
                    9.07005884e-10, -9.04424499e-13, -1.43440860e4, 3.50840928e0],
            },
            NasaPolynomial {
                t_low: 1000.0, t_high: 6000.0,
                a: [3.04848583e0, 1.35172818e-3, -4.85794075e-7,
                    7.88536486e-11, -4.69807489e-15, -1.42661171e4, 6.01709790e0],
            },
        ],
    }
}

/// Returns NASA-7 polynomial data for N₂ (molecular nitrogen).
pub fn n2() -> NasaSpecies {
    NasaSpecies {
        name: "N2".into(),
        molar_mass: 28.014,
        polynomials: vec![
            NasaPolynomial {
                t_low: 200.0, t_high: 1000.0,
                a: [3.53100528e0, -1.23660988e-4, -5.02999433e-7,
                    2.43530612e-9, -1.40881235e-12, -1.04697628e3, 2.96747038e0],
            },
            NasaPolynomial {
                t_low: 1000.0, t_high: 6000.0,
                a: [2.95257637e0, 1.39690040e-3, -4.92631603e-7,
                    7.86010195e-11, -4.60755204e-15, -9.23948688e2, 5.87188762e0],
            },
        ],
    }
}
```

### Where to put these files

| New file | Contents |
|---|---|
| `src/core/nasa.rs` | `NasaPolynomial`, `NasaSpecies`, hard-coded species constructors |
| `src/core/thermo.rs` | `ThermodynamicGas` trait, `FrozenGas` (used in V3) |
| `src/core/mod.rs` | Add `pub mod nasa;` and `pub mod thermo;` |

---

## Practical Path for V2

The minimal V2 implementation avoids a full rewrite of the PM function while still improving
accuracy significantly:

1. **Hard-code NASA-7 coefficients** for 3–4 species relevant to LOX/LH₂ or LOX/RP-1 combustion
   products (H₂O, H₂, CO₂, CO, N₂) — the functions above provide this directly.

2. **Implement `NasaSpecies::cp_over_r()` and `enthalpy_j_per_kg()`** as shown — these are the
   only two methods needed by the energy equation.

3. **Implement a `MixtureGas` struct** that holds a `Vec<(NasaSpecies, f64)>` (species, mass
   fraction) and computes mixture Cp as the mass-fraction-weighted sum of species Cp values.

4. **Compute γ_eff** as the average of γ(T) over the temperature range from T_throat to T_exit:
   ```
   γ_eff = (1/N) · Σᵢ γ(T_throat + i·ΔT)
   ```
   This is a single pass over N = 20 temperature points and takes microseconds.

5. **Use γ_eff in the existing MOC solver** as a drop-in replacement for `config.gamma`. No
   changes to the characteristic mesh, expansion fan, or wall-point calculations are required.

The **full variable-γ implementation** requires:
- Replacing the closed-form `prandtl_meyer(M)` with `prandtl_meyer_numerical(M, T₀, T_min, N)`
- Solving T(M) at every mesh point using bisection
- Propagating T through the mesh alongside the Mach field

This is more work (roughly 3× the code of the simple γ_eff path) but follows the same structural
approach to the mesh. The `ThermodynamicGas` trait methods above encapsulate the extra complexity.

---

## Summary

- Real combustion gases have γ(T) ranging from ~1.15 (near chamber) to ~1.35 (near exit) —
  using a constant 1.4 introduces significant error in both exit Mach number and nozzle contour
- **NASA-7 polynomials**: Cp/R = a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴ (7 coefficients per temperature
  range, typically two ranges per species: 200–1000 K and 1000–6000 K)
- γ(T) = Cp(T) / (Cp(T) - R) — derived directly from the polynomial with no integration required
- The Prandtl-Meyer function must be integrated numerically when γ varies; for γ varying less than
  ~10%, a mean effective γ_eff is an acceptable engineering approximation
- **Minimum viable V2**: compute γ_eff as the temperature-averaged γ over the expansion range and
  substitute into the existing constant-γ solver — immediate accuracy improvement with minimal code
  change
- **Full V2**: modify `GasModel` to be temperature-aware via the `ThermodynamicGas` trait; replace
  the closed-form PM formula with the `prandtl_meyer_numerical` method; solve T(M) via bisection at
  each mesh point
- Coefficient data sources: NASA CEA database (bundled with the CEA executable), NIST WebBook
  (https://webbook.nist.gov), Burcat database, or the Python `cantera` package
