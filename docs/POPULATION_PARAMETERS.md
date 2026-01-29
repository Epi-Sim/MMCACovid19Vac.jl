# Population Parameters Documentation

This document describes the parameters available in the `population_params` section of the configuration file. These parameters define the demographic structure, contact patterns, and static mobility backbone of the metapopulation model.

**Note:** For time-varying stochastic mobility configuration, please refer to [MOBILITY_CONFIGURATION.md](MOBILITY_CONFIGURATION.md).

---

## Demographic Structure

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `G_labels` | Array{String} | Labels for the age strata (e.g., `["Young", "Adult", "Elderly"]`). |
| `σ` | Float | Average household size. Used for calculating household density effects. |
| `ξ` | Float | Density factor (Xi). Controls how population density impacts local contact rates. |

---

## Contact Matrices (Social Mixing)

These parameters define how different age groups interact with each other.

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `C` | Matrix{Float} | **Contact Matrix.** A $G \times G$ matrix where $C_{gh}$ represents the probability or intensity of contact between group $g$ and group $h$. |
| `kᵍ` | Array{Float} | **Average Daily Contacts.** Total expected contacts per individual in each group per day. |
| `kᵍ_h` | Array{Float} | **Home Contacts.** Average contacts occurring within the household. |
| `kᵍ_w` | Array{Float} | **Work/School Contacts.** Average contacts occurring at work or school locations. |

---

## Static Mobility Backbone

These parameters define the baseline movement patterns before any stochastic noise is applied.

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `pᵍ` | Array{Float} | **Mobility Fraction.** The fraction of people in each age group $g$ who travel (leave their home patch) on a typical day. $0.0$ = no travel, $1.0$ = everyone travels. |

*Note:* The actual routing of travelers (where they go) is defined by the `mobility_matrix_filename` (CSV) specified in the `data` section of the config, not in `population_params`. This CSV loads the $R_{ij}$ matrix.

---

## Usage in Configuration

```json
"population_params": {
    "G_labels": ["Young", "Adult", "Elderly"],
    "pᵍ": [0.0, 1.0, 0.2],
    "kᵍ": [12.0, 13.5, 7.0],
    "C": [
        [0.8, 0.1, 0.1],
        [0.1, 0.8, 0.1],
        [0.1, 0.1, 0.8]
    ],
    ...
}
```

```