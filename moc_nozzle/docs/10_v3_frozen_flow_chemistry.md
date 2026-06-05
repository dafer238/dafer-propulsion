# V3 — Frozen-Flow Chemistry and NASA CEA Integration

## Combustion Products Are Not Air

The baseline solver uses `γ = 1.4` (air). Real rocket engines burn propellant to produce hot,
high-pressure combustion products whose chemical composition is nothing like air.

For **LOX/LH₂** at O/F ≈ 6, the combustion products at 3500 K and 10 MPa are roughly:

| Species | Mole fraction |
|---------|--------------|
| H₂O     | 0.55         |
| H₂      | 0.22         |
| OH      | 0.10         |
| O       | 0.04         |
| H       | 0.05         |
| O₂      | 0.04         |

This mixture has γ ≈ 1.14 at chamber conditions — dramatically different from 1.4. Using 1.4
overpredicts exit velocity by ~15%.

The core problem is that triatomic molecules like H₂O have many active internal energy modes at
high temperature: three translational, three rotational, and three vibrational modes (compared to
a diatomic like N₂ which has only two rotational and one vibrational mode). More internal energy
modes means more ways the molecule can store thermal energy, which means a lower fraction of that
energy converts to directed kinetic energy (thrust) during expansion — hence the lower γ.

---

## Frozen Flow vs. Equilibrium Flow

Two limiting thermochemical models bound the real behavior:

**Equilibrium flow**: as the gas expands and cools, chemical reactions proceed to thermodynamic
completion at each point in the nozzle. Species fractions change continuously — for example,
recombination reactions like H + H → H₂ release energy as the gas cools, supplementing the
kinetic energy of the flow. This gives the **maximum theoretical Isp** but requires solving a
nonlinear system of equilibrium equations at every mesh point (expensive and numerically
sensitive).

**Frozen flow**: chemical composition is fixed ("frozen") at some reference condition — usually the
throat. As the gas expands downstream of the throat, the species fractions stay constant, but
temperature-dependent properties (Cp, γ) still vary with T through the NASA polynomial for the
frozen mixture. This is more conservative than equilibrium (lower predicted Isp) and far simpler
to implement.

**Which one to use:**
- Frozen flow at the throat is the standard engineering approximation for nozzle contour design
- The true Isp lies between the frozen and equilibrium values
- For H₂/O₂, the difference is typically 20–40 s in Isp (out of ~450 s total)
- Frozen flow from the throat is the recommended model for V3

The physical justification for freezing at the throat is that the gas expansion time scale in the
diverging nozzle is shorter than the chemical relaxation time scale once temperature drops below
~2000 K. At chamber conditions, reactions are fast and equilibrium holds. At the throat
(M = 1, T ≈ 3000–3500 K), the expansion becomes supersonic and the residence time drops sharply.
Most major recombination reactions effectively stop here.

---

## NASA CEA: What It Is and What It Outputs

NASA CEA (Chemical Equilibrium with Applications) is the industry-standard tool for computing
combustion properties. It was originally developed at NASA Glenn Research Center and is distributed
free of charge. CEA:

1. Takes propellant inputs (fuel type, oxidizer type, oxidizer-to-fuel mass ratio O/F, chamber
   pressure Pc)
2. Solves the equilibrium chemistry at the chamber and optionally at specified expansion ratios
3. Outputs species fractions, T, P, γ, Cp, and mean molar mass M̄ at each condition

**Relevant CEA outputs for the MOC solver:**

| Output | Symbol | Units | Used for |
|--------|--------|-------|----------|
| Combustion temperature | Tc | K | Stagnation temperature input |
| Throat temperature | Tt | K | Throat boundary condition |
| Effective γ at throat | γt | — | Frozen-flow PM function |
| Mean molar mass | M̄ | g/mol | R_specific = R_univ / M̄ |
| Cp at throat | Cp_t | J/(kg·K) | Verification, γ_eff calculation |
| Species mole fractions at throat | xᵢ | — | Input to NASA polynomial mixture |

**How to run CEA:**
- Download from NASA Glenn: https://www.grc.nasa.gov/WWW/CEAWeb/
- It is a Fortran executable that reads a plain-text input file (`.inp`) and writes a plain-text
  output file (`.out`)
- Example minimal input (LOX/LH₂, Pc = 100 bar = 10 MPa, O/F = 6, frozen expansion from throat):

```
problem  case = "LOX-LH2-example"
  rocket frozen  nfz=1
  p(bar) = 100
  sup = 10  ! area expansion ratio Ae/At
  phi = 1.0  ! equivalence ratio (phi=1 → stoichiometric; for H2/O2 this is O/F≈8)
reactants
  fuel H2(L) wt=1.0
  oxid O2(L)  wt=1.0
  o/f = 6.0
output
  massf short transport
end
```

The `frozen nfz=1` keyword tells CEA to freeze composition at the throat (station 1 of the nozzle).
The `sup = 10` keyword requests the solution at an expansion ratio of 10.
The `massf` output keyword requests mass fractions (as opposed to mole fractions).

**Alternatives to running CEA directly:**

| Tool | Language | Notes |
|------|----------|-------|
| `rocketcea` | Python | pip-installable wrapper around CEA; excellent for parameter sweeps; returns structured Python objects |
| `thermo` / `cantera` | Python | Full chemical equilibrium; more flexible but heavier dependency |
| Pre-computed tables | Any | Run CEA over a grid of O/F and Pc; store results as JSON; read in Rust at runtime |
| `pyCEA` | Python | Lightweight pure-Python reimplementation (less accurate, good for prototyping) |

For V3 integration in this Rust project, the **recommended path** is:

1. Run CEA or `rocketcea` in Python once per propellant combination
2. Extract the key output parameters (Tc, Tt, γ_t, M̄, species fractions)
3. Store them in a TOML or JSON config file in the project
4. Read the config into a `CeaData` struct at runtime

This keeps the Rust binary dependency-free while leveraging the validated NASA tool for the
thermochemistry.

---

## The Frozen-Flow Isentropic Relations

For frozen flow with fixed composition but temperature-dependent Cp(T) via NASA polynomials, the
isentropic relations differ from the constant-γ case.

The flow is isentropic along streamlines, so entropy is constant:

```
ds = Cp_mix(T) dT/T - R_mix dP/P = 0
```

Integrating from stagnation state (T₀, P₀) to local state (T, P):

```
∫_{T₀}^{T} Cp_mix(T') / T' dT' = R_mix · ln(P/P₀)
```

This gives P(T) for the frozen mixture during isentropic expansion. Note that the standard formula
T/T₀ = (P/P₀)^{(γ-1)/γ} is only the constant-γ special case obtained when Cp is independent of T.

The area-Mach relation is also modified. In terms of the throat (*) and local conditions:

```
A/A* = (ρ* a*) / (ρ a)   ×   (1/M)
```

where throat conditions are computed from the isentropic integral with M = 1. For a simplified
frozen-flow implementation (acceptable as the V3 starting point):

- Compute γ_eff as the NASA-polynomial-based average γ over the temperature range T_throat to
  T_exit (as described in doc 09)
- Use γ_eff in the standard isentropic area-Mach formula
- Use NASA polynomial h(T) to verify energy conservation in post-processing

The error introduced by using γ_eff instead of the full variable-γ integral is typically less than
1% for LOX/LH₂ and less than 2% for LOX/RP-1 — well within the accuracy of the frozen-flow
assumption itself.

---

## Rust Integration Strategy

### Step 1: Define a `CeaData` configuration struct

This struct is the bridge between the external NASA CEA tool and the Rust MOC solver. It carries
all the combustion information the solver needs without requiring CEA to be linked or called at
runtime. Place it in `src/core/cea.rs`:

```rust
/// Combustion properties extracted from a NASA CEA run.
///
/// Populate this struct from a CEA output file, a rocketcea Python call,
/// or from the pre-computed reference table below. Feed it into `NozzleConfig`
/// to replace the constant-γ assumption with physically correct thermochemistry.
#[derive(Clone, Debug)]
pub struct CeaData {
    /// Human-readable propellant description, e.g. "LOX/LH2 O/F=6.0 Pc=10MPa"
    pub propellant:     String,

    /// Chamber (combustion) temperature [K].
    /// This is the stagnation temperature T₀ for the entire nozzle calculation.
    pub t_chamber:      f64,

    /// Chamber pressure [Pa].
    /// Used to normalize pressure ratios in post-processing.
    pub p_chamber:      f64,

    /// Throat temperature after isentropic expansion to M=1 [K].
    /// Used as the starting temperature for the frozen-flow expansion.
    pub t_throat:       f64,

    /// Effective γ at the throat, as output by CEA.
    /// For frozen flow, this is the most important single number: it sets
    /// the throat boundary condition for the MOC expansion fan.
    pub gamma_throat:   f64,

    /// Estimated γ at the exit plane (useful for post-processing Isp).
    /// May be computed from NASA polynomials at T_exit, or approximated as
    /// gamma_throat + 0.05 as a rough correction.
    pub gamma_exit_est: f64,

    /// Mean molar mass of the frozen mixture [g/mol].
    /// Used to compute R_specific = R_universal / M_mol.
    pub molar_mass:     f64,

    /// Specific heat at constant pressure at throat [J/(kg·K)].
    /// For cross-checking: Cp_throat = gamma_throat * R_specific / (gamma_throat - 1).
    pub cp_throat:      f64,

    /// Specific gas constant of the mixture [J/(kg·K)].
    /// Computed as R_universal / M_mol. Stored here for convenience.
    pub r_specific:     f64,
}

impl CeaData {
    /// Construct from CEA output values, computing r_specific automatically.
    pub fn new(
        propellant: impl Into<String>,
        t_chamber: f64,
        p_chamber: f64,
        t_throat: f64,
        gamma_throat: f64,
        gamma_exit_est: f64,
        molar_mass: f64,
        cp_throat: f64,
    ) -> Self {
        CeaData {
            propellant: propellant.into(),
            t_chamber,
            p_chamber,
            t_throat,
            gamma_throat,
            gamma_exit_est,
            molar_mass,
            cp_throat,
            r_specific: 8314.0 / molar_mass,
        }
    }

    /// Speed of sound at the throat [m/s].
    pub fn throat_sound_speed(&self) -> f64 {
        (self.gamma_throat * self.r_specific * self.t_throat).sqrt()
    }

    /// Characteristic velocity c* [m/s].
    ///
    /// c* = sqrt(R T_chamber / gamma) * (gamma+1)/2)^((gamma+1)/(2*(gamma-1)))
    /// This is the standard formula from Sutton & Biblarz.
    pub fn c_star(&self) -> f64 {
        let g = self.gamma_throat;
        let coeff = ((g + 1.0) / 2.0).powf((g + 1.0) / (2.0 * (g - 1.0)));
        (self.r_specific * self.t_chamber / g).sqrt() / coeff
    }
}
```

### Step 2: Pre-computed reference values for common propellants

The following table provides typical CEA results at Pc = 10 MPa. Use these values to get the
solver running immediately without needing to install and run CEA. For final design work, replace
with your specific O/F and Pc conditions.

| Propellant  | O/F | Tc [K] | Tt [K] | γ_t  | M̄ [g/mol] | Isp_vac [s] |
|-------------|-----|--------|--------|------|-----------|------------|
| LOX/LH₂    | 6.0 | 3540   | 3070   | 1.13 | 12.0      | ~455       |
| LOX/RP-1   | 2.7 | 3670   | 3170   | 1.16 | 21.5      | ~363       |
| LOX/CH₄    | 3.4 | 3550   | 3060   | 1.16 | 19.8      | ~380       |
| N₂O₄/UDMH | 2.6 | 3200   | 2750   | 1.19 | 21.4      | ~340       |

> **Note:** Values are approximate averages from literature. Run CEA for your specific conditions;
> γ and Isp are sensitive to O/F ratio, which has an optimum near maximum temperature.

Rust constructor functions for these reference cases:

```rust
impl CeaData {
    /// LOX/LH₂ at O/F = 6.0, Pc = 10 MPa. Approximate reference values.
    pub fn lox_lh2_of6() -> Self {
        Self::new("LOX/LH2 O/F=6.0 Pc=10MPa",
            3540.0, 10.0e6, 3070.0, 1.13, 1.20, 12.0,
            /*cp_throat*/ 1.13 * (8314.0 / 12.0) / (1.13 - 1.0))
    }

    /// LOX/RP-1 at O/F = 2.7, Pc = 10 MPa. Approximate reference values.
    pub fn lox_rp1_of27() -> Self {
        Self::new("LOX/RP-1 O/F=2.7 Pc=10MPa",
            3670.0, 10.0e6, 3170.0, 1.16, 1.22, 21.5,
            1.16 * (8314.0 / 21.5) / (1.16 - 1.0))
    }

    /// LOX/CH₄ at O/F = 3.4, Pc = 10 MPa. Approximate reference values.
    pub fn lox_ch4_of34() -> Self {
        Self::new("LOX/CH4 O/F=3.4 Pc=10MPa",
            3550.0, 10.0e6, 3060.0, 1.16, 1.22, 19.8,
            1.16 * (8314.0 / 19.8) / (1.16 - 1.0))
    }
}
```

### Step 3: Implement `FrozenGas` that uses `CeaData` + `NasaPolynomial`

`FrozenGas` is the central V3 type. It pairs the CEA macro-level data with NASA polynomial species
data for temperature-dependent property evaluation. Place it in `src/core/thermo.rs`:

```rust
use crate::core::gas::GasModel;
use crate::core::nasa::NasaSpecies;
use crate::core::cea::CeaData;

/// A frozen-flow gas model derived from CEA output.
///
/// Composition is fixed at throat conditions (species list and mass fractions
/// from CEA). NASA polynomials for each species provide Cp(T) and h(T),
/// allowing γ(T) to vary continuously during the expansion.
///
/// For a minimal implementation, only `cea` is required. The `species` field
/// is optional: if empty, the struct falls back to using `cea.gamma_throat`
/// as a constant throughout the expansion.
pub struct FrozenGas {
    pub cea:     CeaData,
    /// Pairs of (species data, mass fraction at throat).
    /// Mass fractions must sum to 1.0.
    pub species: Vec<(NasaSpecies, f64)>,
}

impl FrozenGas {
    /// Mixture Cp(T) in J/(kg·K), weighted by mass fractions.
    ///
    /// For each species i with mass fraction mfᵢ:
    ///   Cp_mix = Σᵢ mfᵢ · Cp_i(T)
    /// where Cp_i(T) = (cp_over_r)_i × R_universal / M_mol_i
    pub fn cp_mix(&self, t: f64) -> f64 {
        self.species.iter().map(|(sp, mf)| {
            let r_sp_i = 8314.0 / sp.molar_mass;  // J/(kg·K) for species i
            sp.cp_over_r(t) * r_sp_i * mf
        }).sum()
    }

    /// Mixture enthalpy h_mix(T) in J/kg, weighted by mass fractions.
    pub fn enthalpy_mix(&self, t: f64) -> f64 {
        self.species.iter().map(|(sp, mf)| {
            sp.enthalpy_j_per_kg(t) * mf
        }).sum()
    }

    /// Local mixture γ(T) = Cp_mix(T) / (Cp_mix(T) - R_mix).
    ///
    /// Falls back to `cea.gamma_throat` if no species data is provided.
    pub fn gamma_at_t(&self, t: f64) -> f64 {
        if self.species.is_empty() {
            return self.cea.gamma_throat;
        }
        let cp = self.cp_mix(t);
        let r  = self.cea.r_specific;
        cp / (cp - r)
    }

    /// Mean effective γ over the temperature range [t_low, t_high].
    ///
    /// Evaluated at `n_points` evenly spaced temperatures. Use this to
    /// get a single γ_eff for the constant-γ approximation in V3.
    pub fn gamma_eff(&self, t_low: f64, t_high: f64, n_points: usize) -> f64 {
        let n = n_points.max(2);
        let dt = (t_high - t_low) / (n - 1) as f64;
        let sum: f64 = (0..n)
            .map(|i| self.gamma_at_t(t_low + i as f64 * dt))
            .sum();
        sum / n as f64
    }

    /// Speed of sound at temperature T [m/s].
    pub fn sound_speed(&self, t: f64) -> f64 {
        (self.gamma_at_t(t) * self.cea.r_specific * t).sqrt()
    }
}

impl GasModel for FrozenGas {
    /// Return γ at the throat as the representative constant for the baseline MOC.
    ///
    /// This is the "minimum viable" V3: simply using the correct γ for the
    /// propellant combination instead of the air value γ = 1.4.
    fn gamma(&self) -> f64 {
        self.cea.gamma_throat
    }

    fn prandtl_meyer(&self, m: f64) -> f64 {
        // Uses constant gamma_throat as the approximation.
        // For full variable-γ, replace with prandtl_meyer_numerical() from
        // ThermodynamicGas (doc 09), passing self.cea.t_throat as T₀.
        let g = self.gamma();
        let a = ((g + 1.0) / (g - 1.0)).sqrt();
        let b = ((g - 1.0) / (g + 1.0) * (m * m - 1.0)).sqrt();
        a * b.atan() - (m * m - 1.0_f64).sqrt().atan()
    }

    fn inverse_prandtl_meyer(&self, nu: f64) -> f64 {
        use crate::utils::root::bisection;
        if nu <= 0.0 { return 1.0; }
        bisection(|m| self.prandtl_meyer(m) - nu, 1.0 + 1e-9, 100.0)
    }
}
```

### Step 4: Feed `CeaData` into `NozzleConfig`

Extend the existing config struct to accept optional CEA data. When present, it overrides the
fallback constant γ:

```rust
/// Configuration for the MOC nozzle solver.
pub struct NozzleConfig {
    /// Fallback constant γ. Used when `cea` is None.
    /// Default: 1.4 (air, for development/testing only).
    pub gamma:         f64,

    /// Area expansion ratio Ae/A* (exit area / throat area).
    pub ae_at:         f64,

    /// Number of characteristic lines in the initial expansion fan.
    /// More lines → finer mesh → higher accuracy → slower computation.
    pub n_chars:       usize,

    /// Throat radius [m]. Sets the physical scale of the nozzle.
    pub throat_radius: f64,

    /// True for an axisymmetric (bell) nozzle; false for a 2D planar nozzle.
    pub axisymmetric:  bool,

    /// NASA CEA combustion data. When Some, replaces `gamma` with
    /// `cea.gamma_throat` and enables NASA polynomial property lookups.
    pub cea:           Option<CeaData>,
}

impl NozzleConfig {
    /// Return the effective γ to use in the MOC solver.
    ///
    /// Prefers CEA throat γ when available; falls back to config.gamma.
    pub fn effective_gamma(&self) -> f64 {
        self.cea.as_ref().map(|c| c.gamma_throat).unwrap_or(self.gamma)
    }
}
```

### Step 5: Use effective γ in the solver

The minimum code change to activate V3 in the existing solver is a single-line substitution in
`src/solver/moc.rs` (or wherever `config.gamma` is first used to construct the `GasModel`):

```rust
// Before (V1):
let gas = PerfectGas::new(config.gamma);

// After (V3 minimum viable):
let gas = PerfectGas::new(config.effective_gamma());

// After (V3 full, using FrozenGas):
let gas: Box<dyn GasModel> = if let Some(cea) = &config.cea {
    Box::new(FrozenGas {
        cea: cea.clone(),
        species: build_species_mixture(cea),  // from NASA polynomial data
    })
} else {
    Box::new(PerfectGas::new(config.gamma))
};
```

The `build_species_mixture` function looks up species from the NASA polynomial database and assigns
mass fractions from the CEA output. For a simple V3 without per-point T tracking, this `FrozenGas`
instance uses `gamma_throat` throughout — still a large improvement over γ = 1.4.

---

## Understanding the "Frozen at Throat" Assumption

It is worth unpacking what "frozen at throat" means physically and what its limitations are:

**Why the throat?** In the converging section and at the chamber, the gas is in near-equilibrium:
pressures and temperatures are high, residence times are long, and chemical reaction rates are fast.
The standard CEA solution gives equilibrium composition here. At the throat, the Mach number
reaches 1 and the flow becomes supersonic. In the diverging nozzle, the expansion is rapid: the
temperature drops hundreds of Kelvin per centimeter, and the reaction time scales for species like
H + OH → H₂O become longer than the flow time scale. Composition freezes.

**What stays constant:** Species mole (or mass) fractions. The amounts of H₂O, H₂, CO₂, CO,
OH, etc. are fixed at their throat values.

**What still varies:** Temperature T, pressure P, density ρ, velocity V — all these change during
the isentropic expansion. And because Cp of each species is temperature-dependent (via NASA
polynomials), the mixture Cp(T) and γ(T) also vary, even with fixed composition. This is the V2
improvement applied within the frozen V3 framework.

**Key frozen-flow outputs from CEA needed for the MOC solver:**

| Quantity | Why needed |
|---------|-----------|
| T_throat | Starting temperature for downstream T(M) calculation |
| γ(T_throat) | Throat PM function and boundary condition |
| Species fractions at throat | Build the NASA polynomial mixture for Cp_mix(T) |
| M̄ at throat | R_specific for the entire expansion |
| Cp(T_throat) | Cross-check: should match γ·R_sp/(γ-1) |

---

## The Specific Impulse Connection

Specific impulse is the ultimate metric of rocket engine performance:

```
Isp = v_e / g₀   [seconds]
```

where v_e is the effective exhaust velocity and g₀ = 9.80665 m/s² is standard gravity.

The effective exhaust velocity from the MOC exit conditions is:

```
v_e = M_exit · a_exit = M_exit · sqrt(γ(T_exit) · R_mix · T_exit)
```

Or equivalently from the energy equation:

```
v_e = sqrt(2 · (h_throat - h_exit))   [frozen flow, energy balance]
```

**Why γ matters for Isp:** Using γ = 1.4 (air) instead of γ ≈ 1.14 (LOX/LH₂) changes the
speed of sound:

```
a_exit = sqrt(γ · R · T_exit)
```

The ratio of sound speeds is:

```
a_exit(γ=1.4) / a_exit(γ=1.14) = sqrt(1.4 / 1.14) ≈ 1.11
```

So using the wrong γ causes an **11% error in Isp** from the sound-speed term alone, before even
accounting for the change in exit Mach number caused by the different PM function. In absolute
terms, for a LOX/LH₂ engine with true Isp_vac ≈ 450 s, this 11% error means ~50 s error in Isp —
which is the difference between a competitive engine and an uncompetitive one.

Breaking down the V1 → V3 improvement:

| Model | γ used | Isp_vac [s] (LOX/LH₂, Ae/At=40) | Error vs. CEA |
|-------|--------|----------------------------------|--------------|
| V1: constant air γ | 1.40 | ~500 | +10–12% |
| V2: γ_eff from NASA poly | 1.16 | ~455 | +1–2% |
| V3: frozen flow, γ_t from CEA | 1.13 | ~453 | < 1% |
| Full equilibrium (reference) | variable | ~455 | 0% (reference) |

The V3 frozen-flow value is slightly conservative (lower than equilibrium) because it neglects
recombination energy release in the diverging section.

---

## Summary

- Real rocket propellant combustion products have γ ≈ 1.13–1.20, not 1.4; using the wrong γ
  introduces ~10% error in Isp and significantly distorts the nozzle contour
- **Frozen flow**: composition is fixed ("frozen") at the throat; γ(T) still varies downstream via
  NASA polynomials for the frozen mixture — simpler than equilibrium and accurate to ~1–2% in Isp
- **NASA CEA** provides chamber and throat conditions; it is best run externally (as a Fortran
  executable or via the Python `rocketcea` wrapper), with results stored in a TOML/JSON config file
  that the Rust solver reads at startup
- **`CeaData` struct**: carries Tc, Tt, γ_t, M̄, Cp_t, and R_specific — the complete set of
  combustion properties needed by the MOC solver
- **`FrozenGas` struct**: wraps `CeaData` + NASA polynomial species data; implements `GasModel`
  using `gamma_throat` as the constant approximation, with `gamma_at_t(T)` available for the
  full variable-γ path
- **Minimum viable V3**: replace `gamma = 1.4` with `gamma = cea.gamma_throat` from the CEA
  output — one line of code change, immediate 8–10% improvement in Isp accuracy, no algorithmic
  changes to the MOC solver required
- **Full V3**: use `gamma_at_t(T)` everywhere in the solver; integrate T(M) from the energy
  equation at each mesh point; use numerical PM function (see doc 09); propagate T through the
  characteristic mesh alongside M
- **Isp improvement from V3 alone**: 5–10% over using γ = 1.4, closing most of the gap to a
  full equilibrium calculation at a fraction of the complexity
