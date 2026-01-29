# MMCACovid19Vac

MMCACovid19Vac is a Julia package for simulating epidemic dynamics in a metapopulation, with built-in support for modeling vaccination campaigns. The package is based on two previous implementations:
* MMCAcovid19 Julia package [https://docs.juliahub.com/MMCAcovid19/]([url](https://docs.juliahub.com/MMCAcovid19/))
* MMCA_with_vaccination [https://github.com/PGcastioni/MMCA_with_vaccination/]([url](https://github.com/PGcastioni/MMCA_with_vaccination/))

# Brief history

## MMCAcovid19
This package [MMCAcovid19](https://github.com/jtmatamalas/MMCAcovid19), written in the [Julia](https://julialang.org) language, implements the epidemic model for COVID-19 developed by a group of researchers from [Universitat Rovira i Virgili](https://www.urv.cat) and [Universidad de Zaragoza](http://unizar.es) [[1](#References-1)]. The model makes use of a Microscopic Markov Chain Approach (MMCA) to describe mathematically the dynamics of a so-called metapopulation model of epidemic spreading [[2-4](#References-1)]. 

## MMCA_with_vaccination
The original was later extended by Piergiorgio Castioni to allow modelling of the effect of vaccines, herd immunity and reinfections. Each of the previous compartments has been triplicated to account for the fact that there are vaccinated and unvaccinated individuals. Having received the vaccine changes some of the transition probabilities related to the most negative aspects of the disease, such as transmission, hospitalization and death.

## MMCACovid19Vac

This project is an attempt to extend the functionalities of the previous version. Currently, this version provides a standard configuration format for defining a model and provides a simple set of command line scripts to generate configuration templates and run the model.

## Configuration & Parameters

The model is configured using a JSON file containing several sections. Detailed documentation for each parameter set is available in the `docs/` directory:

*   **[Epidemic Parameters](docs/EPIDEMIC_PARAMETERS.md)**: Disease progression, transmission rates, and clinical outcomes.
*   **[Population Parameters](docs/POPULATION_PARAMETERS.md)**: Demographics, contact matrices, and static mobility.
*   **[Mobility Configuration](docs/MOBILITY_CONFIGURATION.md)**: Time-varying stochastic mobility and noise injection.
*   **[NPI Parameters](docs/NPI_PARAMETERS.md)**: Non-Pharmaceutical Interventions (lockdowns, social distancing).

### Quick Start

This package is designed to be used as a library within your Julia scripts.

#### 1. Running a Simulation

Create a Julia script (e.g., `run_sim.jl`) with the following structure:

```julia
using MMCACovid19Vac
using JSON
using Dates
using CSV
using DataFrames

# 1. Load Configuration
config = JSON.parsefile("test/data/config_test.json") # Adjust path as needed
data_dict = Dict(config["data"])
pop_params_dict = Dict(config["population_params"])
epi_params_dict = Dict(config["epidemic_params"])
npi_params_dict = Dict(config["NPI"])
vac_params_dict = Dict(config["vaccination"])

# 2. Load Data Frames
# Note: Paths should be relative to where you run the script
metapop_df = CSV.read(joinpath("test/data", data_dict["metapopulation_data_filename"]), DataFrame, types=Dict(:id => String))
network_df = CSV.read(joinpath("test/data", data_dict["mobility_matrix_filename"]), DataFrame)

# 3. Initialize Parameters
G = 3 # Number of age groups
M = size(metapop_df, 1) # Number of patches
T = (Date(config["simulation"]["end_date"]) - Date(config["simulation"]["start_date"])).value + 1

# Initialize Structures
population = init_pop_param_struct(G, M, map(String, pop_params_dict["G_labels"]), pop_params_dict, metapop_df, network_df)
epi_params = init_epi_parameters_struct(G, M, T, map(String, pop_params_dict["G_labels"]), epi_params_dict)
npi_params = init_NPI_parameters_struct("test/data", npi_params_dict, get(data_dict, "kappa0_filename", nothing), Date(config["simulation"]["start_date"]))

# 4. Set Initial Conditions
# (Example: 1 exposed individual in patch 1, age group 2)
E₀ = zeros(Float64, G, M)
E₀[2, 1] = 1.0
set_initial_conditions!(epi_params, population, zeros(G,M), E₀, zeros(G,M), zeros(G,M), zeros(G,M), zeros(G,M))

# 5. Run Simulation
# Setup vaccination schedule (example: none)
tᵛs = [Int(config["vaccination"]["start_vacc"]), Int(config["vaccination"]["start_vacc"]) + Int(config["vaccination"]["dur_vacc"]), T]
ϵᵍs = zeros(G, 3) 

println("Starting simulation...")
run_epidemic_spreading_mmca!(epi_params, population, npi_params, tᵛs, ϵᵍs; verbose = true)
println("Simulation complete.")

# 6. Save Results
save_simulation_hdf5(epi_params, population, "output.h5")
```

#### 2. Generating a Configuration Template

You can generate a default JSON configuration file to start customizing your model:

```julia
using MMCACovid19Vac
using JSON

# Create a template for 3 age groups
config = create_config_template(3)

# Save to file
open("my_config.json", "w") do f
    JSON.print(f, config, 4) # Pretty print with indentation
end
```

### Installation

```julia
using Pkg
Pkg.add(url="https://github.com/YourRepo/MMCACovid19Vac.jl")
```


1. Alex Arenas, Wesley Cota, Jesús Gómez-Gardeñes, Sergio Gómez, Clara Granell, Joan T. Matamalas, David Soriano-Paños and Benjamin Steinegger: Modeling the spatiotemporal epidemic spreading of COVID-19 and the impact of mobility and social distancing interventions, _Physical Review X_ **10** (2020) 041055 ([doi](https://doi.org/10.1103/PhysRevX.10.041055))

2. Sergio Gómez, Alex Arenas, Javier Borge-Holthoefer, Sandro Meloni and Yamir Moreno: Discrete-time Markov chain approach to contact-based disease spreading in complex networks, _Europhysics Letters_ **89** (2010) 38009 ([doi](https://doi.org/10.1209/0295-5075/89/38009))

3. Jesús Gómez-Gardeñes, David Soriano-Paños and Alex Arenas: Critical regimes driven by recurrent mobility patterns of reaction-diffusion processes in networks, _Nature Physics_ **14** (2018) 391–395 ([doi](https://doi.org/10.1101/2020.03.21.20040022))

4. David Soriano-Paños, L. Lotero, Alex Arenas and Jesús Gómez-Gardeñes: Spreading processes in multiplex metapopulations containing different mobility networks, _Physical Review X_ **8** (2018) 031039 ([doi](https://doi.org/10.1103/PhysRevX.8.031039))
