# MOC Nozzle Design Procedure

## Minimum-Length Nozzle (MLN)

A minimum-length nozzle (MLN) achieves the design exit Mach number in the shortest possible axial distance. The key insight is:

- The throat lip causes an expansion fan (Prandtl-Meyer) from M=1 to some higher Mach
- This expansion fan is the "initial data" for the MOC
- The nozzle wall is designed so that after all reflections, the exit flow is exactly uniform and parallel at M_e

**Design condition for MLN: θ_max = ν_e / 2**

Where:
- ν_e = ν(M_e) = Prandtl-Meyer angle at exit Mach number
- θ_max = maximum wall angle (at the throat lip)

**Why θ_max = ν_e / 2?**

In the simple-wave expansion at the throat, K⁻ = θ − ν = 0 for all initial characteristics (they all start from the sonic throat where θ=0, ν=0). So ν = θ in the initial fan. At the wall, the flow has turned by θ_max and accelerated to ν = θ_max. After reflecting off the axis and reaching the exit, the total accumulated ν must equal ν_e = θ_max + θ_max = 2·θ_max. Hence:

```
θ_max = ν_e / 2
```

---

## Step 1 — Determine Exit Mach Number

**Given:** `ae_at` (exit-to-throat area ratio), γ

**Compute:**
- M_e using `mach_from_area_ratio(ae_at)` — inversion of the isentropic area-Mach relation, supersonic branch
- ν_e = ν(M_e) using `prandtl_meyer(M_e)` (in radians)
- θ_max = ν_e / 2

**Example** for `ae_at = 10.0`, γ = 1.4:

| Quantity | Value |
|----------|-------|
| M_e | ≈ 3.96 |
| ν_e | ≈ 1.318 rad |
| θ_max | ≈ 0.659 rad ≈ 37.7° |

---

## Step 2 — Generate the Initial Data Line

The initial data line represents the first characteristic line just downstream of the throat. In the simple-wave expansion fan, all n initial nodes satisfy **K⁻ = θ − ν = 0**.

Nodes are ordered **from wall to axis** (index 0 = wall, index n−1 = axis):

For k = 0 to n−1:
- t = k / (n − 1) as f64
- θ_k = (1.0 − t) · θ_max  → θ_max at wall (k=0), 0 at axis (k=n−1)
- ν_k = θ_k  (simple wave: K⁻ = 0)
- M_k = ν⁻¹(ν_k)
- x_k = 0.0
- y_k = (1.0 − t) · r_t  → r_t at wall, 0 at axis

**Special treatment:** for the axis node, M → 1.0 and ν → 0, θ → 0.

Also record: the initial wall node **W₀** at (0, r_t) with θ_W = θ_max. This is the starting wall point for the contour.

**Example table** for n=4, ae_at=10.0, γ=1.4 (approximate):

| k | y | θ (rad) | ν (rad) | M |
|---|---|---------|---------|---|
| 0 (wall) | 1.000 | 0.659 | 0.659 | ~1.95 |
| 1 | 0.667 | 0.439 | 0.439 | ~1.57 |
| 2 | 0.333 | 0.220 | 0.220 | ~1.24 |
| 3 (axis) | 0.000 | 0.000 | 0.000 | 1.00 |

---

## Step 3 — March the Characteristic Mesh

The mesh marching loop runs for n rows. In each row:

### 3a. Interior nodes (between axis and wall)

For each adjacent pair (k, k+1) in the current row where k=0 is the wall-side node:
- Node at index k is the **upper/wall-side** node → provides K⁺_R = θ_k + ν_k
- Node at index k+1 is the **lower/axis-side** node → provides K⁻_L = θ_{k+1} − ν_{k+1}
- Compute the new interior node using the interior-point formula:

```
θ_P = (K⁺_R + K⁻_L) / 2
ν_P = (K⁺_R − K⁻_L) / 2
M_P = ν⁻¹(|ν_P|)
```

Position at the intersection of the C⁻ characteristic from the upper node and the C⁺ characteristic from the lower node:

```
s1 = tan(θ_lower + μ_lower)   # C⁺ slope from lower node
s2 = tan(θ_upper − μ_upper)   # C⁻ slope from upper node
x_P = (y_upper − y_lower + s1·x_lower − s2·x_upper) / (s1 − s2)
y_P = y_lower + s1·(x_P − x_lower)
```

### 3b. Axis node (from the last = axis-side node in current row)

The C⁻ characteristic from the axis-side node J reaches the axis. By symmetry, θ = 0 at the axis:

```
K⁺ preserved along C⁻:   K⁺_J = θ_J + ν_J
At axis:   θ = 0  →  ν_axis = K⁺_J
M_axis = ν⁻¹(ν_axis)
```

Position along the C⁻ characteristic from J to y = 0:

```
slope = tan(θ_J − μ_J)
x_axis = x_J − y_J / slope
y_axis = 0.0
```

### 3c. Wall node (from the first = wall-side node J and previous wall node W_prev)

The wall is a streamline. The new wall state comes from the intersection of:
- The C⁻ characteristic through the adjacent interior node J
- The C⁺ characteristic propagating along the wall from W_prev

```
K⁺_J    = θ_J + ν_J             # from C⁻ through J
K⁻_prev = θ_W_prev − ν_W_prev   # from C⁺ along wall

θ_W = (K⁺_J + K⁻_prev) / 2
ν_W = (K⁺_J − K⁻_prev) / 2
```

Position at the intersection of the C⁻ from J and the C⁺ from W_prev:

```
s_J = tan(θ_J − μ_J)                   # C⁻ slope from j_node
s_W = tan(θ_W_prev + μ_W_prev)         # C⁺ slope from w_prev
x_W = (y_J − y_W_prev + s_W·x_W_prev − s_J·x_J) / (s_W − s_J)
y_W = y_J + s_J·(x_W − x_J)
```

### New row assembly

Collect nodes in wall-to-axis order:

```
[wall_node, interior_nodes..., axis_node]
```

---

## Step 4 — Exit Condition

The mesh naturally terminates after n rows (for an MLN with n initial characteristics). At this point:

- The axis node should have **ν ≈ ν_e** (the design PM angle)
- **θ ≈ 0** (parallel flow at exit)
- The last wall node should also have θ → 0

Add a validation check:

```
assert |ν_final_axis − ν_e| < tolerance   (e.g., 1e-3 rad)
```

If this assertion fails, it indicates a bug in the characteristic marching or the initial data generation.

---

## Step 5 — Extract the Wall Contour

The nozzle divergent wall contour is the sequence of wall nodes accumulated during the march:

```
W₀ (throat lip, x=0, y=r_t) → W₁ → W₂ → ... → W_n (exit wall point)
```

These define the inner surface of the diverging nozzle section from throat to exit.

In `geometry/wall.rs`, `extract_wall()` must be rewritten to simply return these accumulated wall nodes in order (see the implementation guide for details). The current heuristic of filtering `y >= 0.5` is incorrect and must be removed.

At the exit:
- The nozzle height = y of the last wall node
- For **2D planar** geometry: A_e/A* ≈ y_exit / r_t
- For **axisymmetric** geometry (using radii): A_e/A* ≈ (y_exit / r_t)²

For the unit throat radius r_t = 1.0, the 2D planar exit area ratio is simply y_exit.

---

## Mermaid Flowchart — Design Procedure

```mermaid
flowchart TD
    A[Input: gamma, ae_at, n_chars] --> B[Compute M_e from ae_at]
    B --> C[Compute nu_e = PM\(M_e\)]
    C --> D[theta_max = nu_e / 2]
    D --> E[Generate initial data line\nn nodes, K- = 0]
    E --> F{For each row\nrow = 1 to n}
    F --> G[Compute interior nodes]
    F --> H[Compute axis node\ntheta=0, nu=K+_J]
    F --> I[Compute wall node\nstreamline BC]
    G & H & I --> J[New characteristic row]
    J --> F
    F -- done --> K[Validate exit condition\nnu_axis approx nu_e]
    K --> L[Extract wall contour\nfrom wall nodes]
    L --> M[Output: wall contour points]
```

---

## Summary

- **MLN design condition:** θ_max = ν_e / 2
- **Initial data:** n nodes with θ = ν (simple wave, K⁻ = 0), ordered wall-to-axis
- **Node ordering:** wall-to-axis in the array — index 0 is wall side, index n−1 is axis side
- **Each row:** interior points (using adjacent pairs) + axis BC (θ=0 symmetry) + wall BC (streamline condition)
- **Exit:** axis node reaches ν ≈ ν_e after n rows; θ → 0 everywhere at the exit plane
- **Wall contour** = accumulated wall nodes from all rows, starting at the throat lip W₀
