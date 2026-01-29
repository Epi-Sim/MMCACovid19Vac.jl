# Mobility Configuration Guide

This document details how to configure time-varying mobility in the MMCACovid19Vac model.

## Overview

In the standard MMCA model, the mobility matrix $R_{ij}$ (fraction of people living in $i$ who visit $j$) is static. This extension allows $R_{ij}$ to vary at each timestep $t$, simulating the stochastic nature of daily human movement while maintaining physical consistency (population conservation).

## Configuration Parameters

These parameters are set in the `population_params` section of your JSON config file.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `mobility_variation_type` | String | `"none"` | The method for generating mobility variation. |
| `mobility_sigma_O` | Float | `0.0` | Noise level for total outflows (Origin strength). |
| `mobility_sigma_D` | Float | `0.0` | Noise level for total inflows (Destination attraction). |
| `mobility_rng_seed` | Int | `null` | Seed for reproducibility. |

### detailed Description

1.  **`mobility_variation_type`**:
    *   `"none"`: Mobility is static. The model uses the input `R_mobility_matrix.csv` for the entire simulation.
    *   `"ipfp_simple"`: Enables the stochastic engine. At each timestep $t$, the model generates a perturbed matrix $R(t)$ based on the baseline matrix $R_{base}$.

2.  **`mobility_sigma_O` & `mobility_sigma_D`**:
    *   These control the variance of the lognormal noise applied to the row sums (origins) and column sums (destinations) of the mobility flux matrix.
    *   Higher values (e.g., `0.2` vs `0.05`) result in more erratic daily changes in how many people leave a patch or visit a specific destination.
    *   The noise is unbiased (mean 1.0) so the *average* mobility over time remains close to the baseline.

3.  **`mobility_rng_seed`**:
    *   If provided, the sequence of random matrices will be identical across runs, allowing for reproducible debugging or sensitivity analysis.

## How it Works (Algorithm)

When `mobility_variation_type` is set to `"ipfp_simple"`:

1.  **Perturbation**: At each step $t$, we take the baseline flows $F_{ij} = N_i R_{ij}$ and calculate the baseline marginals (total outflows $O_i$ and total inflows $D_j$).
2.  **Noise Injection**: We perturb these marginals:
    $$ O_i(t) = O_i \cdot \text{Lognormal}(0, \sigma_O) $$
    $$ D_j(t) = D_j \cdot \text{Lognormal}(0, \sigma_D) $$
3.  **Consistency**: We adjust $D(t)$ slightly so that $\sum O(t) = \sum D(t)$ (global population conservation).
4.  **Reconstruction (IPFP)**: We use the **Iterative Proportional Fitting Procedure (RAS Algorithm)** to find a new matrix $R_{ij}(t)$ that:
    *   Matches the new noisy marginals $O(t)$ and $D(t)$.
    *   Preserves the "structure" (zero/non-zero pattern and relative weights) of the baseline matrix.
5.  **Update**: The model recalculates all density-dependent factors (effective population $\tilde{n}_i$, normalization $z^g$) using this new $R_{ij}(t)$ before computing infection probabilities for the day.

## Limitations

*   **Stochastic Only**: Currently, you cannot provide a CSV or file-based sequence of pre-computed mobility matrices (e.g., `R_t1.csv`, `R_t2.csv`). The variation is generated mathematically at runtime.
*   **Computational Cost**: Enabling this feature adds a small computational overhead per timestep due to the IPFP matrix balancing step.
