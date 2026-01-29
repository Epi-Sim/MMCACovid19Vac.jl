# Epidemic Parameters Documentation

This document describes the parameters available in the `epidemic_params` section of the configuration file. These parameters control the disease progression, transmission rates, and clinical outcomes for different age strata.

---

## Core Transmission Parameters

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `βᴵ` | Float | Infectivity (transmission probability) of symptomatic infected individuals. |
| `βᴬ` | Float | Infectivity of asymptomatic infected individuals. Often scaled relative to `βᴵ` (e.g., `βᴬ = 0.5 * βᴵ`). |
| `scale_β` | Float | **Optional helper.** If `βᴬ` is not explicitly set, it is calculated as `scale_β * βᴵ`. |

---

## Disease Progression Rates (Transition Probabilities)

These parameters define the probability of moving between compartments (e.g., Exposed -> Infected) at each timestep. Values are generally `1 / (average duration in days)`.

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `ηᵍ` | Array{Float} | **Exposed rate.** Probability of becoming infectious (E -> A) per day. |
| `αᵍ` | Array{Float} | **Asymptomatic infectious rate.** Probability of becoming symptomatic (A -> I) per day. |
| `μᵍ` | Array{Float} | **Infectious rate.** Probability of recovery or hospitalization (I -> R or I -> PH/PD) per day. |

---

## Clinical Outcome Probabilities (Stratified)

These parameters are arrays of size $G$ (number of age groups), defining the severity of the disease for each group.

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `θᵍ` | Array{Float} | **Direct Death Probability.** Probability that an infected individual dies without hospitalization (I -> PD). |
| `γᵍ` | Array{Float} | **ICU/Hospitalization Probability.** Probability that an infected individual requires hospitalization (I -> PH). |
| `ωᵍ` | Array{Float} | **Fatality in ICU.** Probability that a hospitalized patient dies (PH -> HD). |

---

## Hospitalization Dynamics

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `ζᵍ` | Array{Float} | **Pre-deceased rate.** Rate of progression to death for non-hospitalized cases (PD -> D). |
| `λᵍ` | Array{Float} | **Pre-hospitalized rate.** Rate of progression from pre-hospital to ICU admission (PH -> HR/HD). |
| `ψᵍ` | Array{Float} | **Death rate in ICU.** Rate of death for fatal ICU cases (HD -> D). |
| `χᵍ` | Array{Float} | **ICU Discharge rate.** Rate of recovery for non-fatal ICU cases (HR -> R). |

---

## Immunity and Reinfection

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `Λ` | Float | **Waning Vaccine Immunity.** Probability of losing vaccine-acquired immunity per day (V -> S). |
| `Γ` | Float | **Waning Natural Immunity.** Probability of losing disease-acquired immunity per day (R -> S). |

---

## Vaccination Efficacy

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `rᵥ` | Array{Float} | **Infection Prevention.** Relative risk reduction of *acquiring* infection for each vaccination status (0=no effect, 1=perfect immunity). |
| `kᵥ` | Array{Float} | **Transmission Prevention.** Relative risk reduction of *transmitting* infection if infected (0=no effect, 1=non-infectious). |

---

## Risk Reduction Factors (Configuration Helpers)

The configuration loader often uses these multipliers to derive the final stratified probabilities for vaccinated groups or specific scenarios.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `risk_reduction_dd` | Float | `0.0` | Multiplier for direct death probability (`θᵍ`) for vaccinated/protected groups. |
| `risk_reduction_h` | Float | `0.1` | Multiplier for hospitalization probability (`γᵍ`). |
| `risk_reduction_d` | Float | `0.05` | Multiplier for ICU fatality probability (`ωᵍ`). |

---

## Usage in Configuration

```json
"epidemic_params": {
    "βᴵ": 0.09,
    "scale_β": 0.5,
    "ηᵍ": [0.275, 0.275, 0.275],
    "αᵍ": [0.65, 0.65, 0.65],
    "μᵍ": [0.3, 0.3, 0.3],
    "θᵍ": [0.0, 0.0, 0.001],
    "γᵍ": [0.01, 0.05, 0.15],
    ...
}
```
