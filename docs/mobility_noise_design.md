# Mobility Variation Design (IPFP/RAS)

## Overview

This document outlines the design for introducing time-varying, structured stochasticity into the mobility network using **Iterative Proportional Fitting Procedure (IPFP)**.

This design strictly adheres to the **Arenas et al. (Phys. Rev. X, 2020)** framework, ensuring that perturbations to the routing matrix $R_{ij}(t)$ are physically consistent with the model's equations for mixing and density.

## Core Concept

The MMCA model drives epidemics through two coupled mechanisms:
1.  **Routing (Eq. 12):** $\Pi_i$ mixes home risk with destination risk based on $R_{ij}$.
2.  **Mixing (Eq. 13-18):** $P_j$ (risk at destination) depends on local effective density $\tilde{n}_j$ and incoming flows $n_{j \to i}$.

**Crucial Constraint:** Because $R_{ij}$ appears in both routing and mixing equations, we cannot simply perturb it in the routing step. We must **re-calculate all population auxiliaries** ($\tilde{n}$, mobility flows, normalization factors) at every timestep where $R_{ij}$ changes.

## 1. Data Structures

We will maintain the **Baseline** (static) mobility and generate a **Transient** (noisy) mobility for the current timestep.

### New Parameters (`Population_Params`)

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `mobility_variation_type` | `String` | `"none"` (default), `"ipfp_simple"` |
| `mobility_sigma_O` | `Float64` | Noise magnitude for Origin totals (outflows). |
| `mobility_sigma_D` | `Float64` | Noise magnitude for Destination totals (inflows). |
| `mobility_rng_seed` | `Union{Int, Nothing}` | Seed for reproducibility. |

*Note: `Rᵢⱼ` in `Population_Params` will serve as the **immutable baseline** ($B_{ij}$). We will not overwrite it.*

## 2. Algorithms

### A. IPFP Solver (Sparse)

We use IPFP to generate a valid stochastic matrix $R(t)$ from the baseline $B$ and noisy marginal targets $O(t), D(t)$.

```julia
function ipfp_sparse(edgelist, B, O_target, D_target; max_iter=20)
    F = copy(B) # Start with baseline structure
    # ... Iterative row/col scaling ...
    # Final step: Enforce strict row-stochasticity (sum_j R_ij = 1.0)
    return F_normalized
end
```

### B. The "Ghost Population" Fix (Auxiliary Updates)

The standard `update_prob!` only updates population params when a confinement event ($t=t_c$) occurs. We must inject a new update routine that runs **every timestep**.

**Standard MMCA Flow:**
1.  `p^g_eff` determined by confinement.
2.  `R_ij` (static).
3.  `update_population_params!`: Computes $n^{eff}, mobility^g, z^g, norm^g$.
4.  `compute_P!`: Uses auxiliaries to get $P$.
5.  `update_prob!`: Uses $R$ and $P$ to get $\Pi$.

**Noisy MMCA Flow:**
1.  **Generate $R(t)$**: Apply IPFP to baseline `population.Rᵢⱼ`.
2.  **Update Auxiliaries**: Call `compute_effective_population!` using $R(t)$ and current $p^g_{eff}$.
    *   *Note: We must be careful not to permanently overwrite the static `population` fields if we want to revert to baseline, or we simply accept `population` represents the "current state". Given the loop structure, updating `population` in-place is the most efficient and standard way, as long as we use the immutable `R_base` for the next step's generation.*
3.  **Compute $P$ & $\Pi$**: Proceed as normal.

## 3. Implementation Plan

### Step 1: Secure the Baseline
In `Population_Params`, the `Rᵢⱼ` field is currently the only storage.
*   **Action:** Add `Rᵢⱼ_base::Array{Float64, 1}` to `Population_Params`.
*   **Init:** In constructors, set `Rᵢⱼ_base = copy(Rᵢⱼ)`.

### Step 2: The `update_mobility!` Function
Create a new function in `src/mmca.jl`:

```julia
function update_mobility!(population::Population_Params, t::Int)
    # 1. Check if noise is enabled
    if population.mobility_variation_type == "none"
        # Optional: Reset to baseline if we were noisy before?
        # Or just do nothing if we assume R is static.
        return
    end

    # 2. Generate Noisy Marginals
    # O(t), D(t) based on population.mobility_sigma...
    
    # 3. Solve IPFP
    R_new = ipfp_sparse(population.edgelist, population.Rᵢⱼ_base, O_t, D_t)
    
    # 4. Update Population State
    population.Rᵢⱼ .= R_new  # Update the "current" R
    
    # 5. Re-calculate Density/Mixing Auxiliaries
    update_population_params!(population) 
end
```

### Step 3: Modify the Main Loop (`update_prob!`)

In `src/mmca.jl`, inside `update_prob!`:

```julia
# ... inside update_prob! ...

# EXISTING: Handle Confinement (Intervention)
if tᶜ == t
    # updates p_eff, etc.
end

# NEW: Handle Stochastic Mobility
# Always run this if noise is enabled, regardless of confinement
if population.mobility_variation_type != "none"
    update_mobility!(population, t)
end

# ... proceed to compute_P! ...
```

## 4. Alternatives Rejected

**Varying `p^g(t)` only:**
*   *Idea:* Use existing `κ₀` mechanism to vary mobility volume.
*   *Verdict:* **Rejected.** `p^g` is a vector of size $G$ (global age strata). It scales mobility uniformly for the entire system. It cannot simulate spatial heterogeneity (e.g., "Region A isolates, Region B travels") or destination changes.

## 5. Performance Considerations
*   **Cost:** IPFP is linear in edges $O(E)$. `update_population_params!` is $O(E \cdot G)$. Total overhead is low for typical metapopulations (e.g., municipal/regional graphs).
*   **Optimization:** If noise is pre-calculated (time-series file) rather than stochastic, we can load $R(t)$ directly. For stochastic noise, on-the-fly IPFP is efficient enough.