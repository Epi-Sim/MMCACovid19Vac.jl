using Test
using DataFrames
using MMCACovid19Vac
using JSON
using Dates
using CSV

## -----------------------------------------------------------------------------
## Population parameters
## -----------------------------------------------------------------------------

G = 3
M = 5
T = 10
V = 3

nᵢᵍ = [ 4995.0   9875.0  14970.0   30010.0   40326.0
       30107.0  59630.0  90009.0  179745.0  239983.0
       15145.0  29827.0  45086.0   90266.0  120026.0]

C = [0.5980 0.3849 0.0171
     0.2440 0.7210 0.0350
     0.1919 0.5705 0.2376]

edgelist = [1  1; 1  2; 1  3; 1  5; 2  1; 2  2; 2  3; 2  4;
            2  5; 3  1; 3  2; 3  3; 3  4; 3  5; 4  1; 4  3;
            4  4; 4  5; 5  1; 5  2; 5  3; 5  5]

Rᵢⱼ = [0.3288; 0.0905; 0.0995; 0.4812; 0.3916; 0.2213; 0.1052; 0.2775;
       0.0044; 0.0233; 0.5205; 0.0117; 0.0807; 0.3638; 0.5156; 0.0579;
       0.0218; 0.4047; 0.3081; 0.2862; 0.0621; 0.3436]

kᵍ = [11.8, 13.3, 6.6]
kᵍ_h = [3.15, 3.17, 3.28]
kᵍ_w = [1.72, 5.18, 0.0]
pᵍ = [0.0, 1.0, 0.05]
sᵢ = [10.6, 23.0, 26.6, 5.7, 61.6]
ξ = 0.01
σ = 2.5

population = Population_Params(G, M, nᵢᵍ, kᵍ, kᵍ_h, kᵍ_w, C, pᵍ, edgelist, Rᵢⱼ, sᵢ, ξ, σ)

## Test population parameters
@testset "Population_Params" begin
    println(pwd())
    @test population.G == G
    @test population.M == M
    @test population.nᵢᵍ == nᵢᵍ
    @test population.kᵍ == kᵍ
    @test population.kᵍ_h == kᵍ_h
    @test population.kᵍ_w == kᵍ_w
    @test population.C == C
    @test population.pᵍ == pᵍ
    @test population.edgelist == edgelist
    @test population.Rᵢⱼ == Rᵢⱼ
    @test population.sᵢ == sᵢ
    @test population.ξ == ξ
    @test population.σ == σ

    @test length(population.nᵢ) == M
    @test length(population.Nᵍ) == G
    @test length(population.nᵢ_eff) == M
    @test length(population.zᵍ) == G
    @test size(population.nᵢᵍ_eff) == (G, M)
    @test size(population.normᵍ) == (G, M)
    @test size(population.mobilityᵍ) == (G, length(Rᵢⱼ))

    @test population.N == Int64(sum(nᵢᵍ))
    @test population.nᵢ ≈ sum(nᵢᵍ, dims=1)[1, :]  atol=0.0001
    @test population.Nᵍ ≈ sum(nᵢᵍ, dims=2)[:, 1]  atol=0.0001

    @test sum(population.nᵢᵍ_eff) ≈ sum(nᵢᵍ)  atol=0.0001
    @test population.nᵢ_eff ≈ sum(population.nᵢᵍ_eff, dims=1)[1, :]  atol=0.0001

end


## -----------------------------------------------------------------------------
## Epidemic parameters
## -----------------------------------------------------------------------------

βᴵ = 0.075
βᴬ = 0.5 * βᴵ
ηᵍ = [1/2.444, 1/2.444, 1/2.444]
αᵍ = [1/5.671, 1/2.756, 1/2.756]
μᵍ = [1/1.0, 1/3.915, 1/3.915]
θᵍ = [0.0, 0.008, 0.047]
γᵍ = [0.0003, 0.003, 0.026]
ζᵍ = [1/7.084, 1/7.084, 1/7.084]
λᵍ = [1/4.084, 1/4.084, 1/4.084]
ωᵍ = [0.3, 0.3, 0.3]
ψᵍ = [1/7.0, 1/7.0, 1/7.0]
χᵍ = [1/20.0, 1/20.0, 1/20.0]

Λ = 0.02
Γ = 0.01
rᵥ = [0.0, 0.6, 0.0]
kᵥ = [0.0, 0.4, 0.0]

risk_reduction_dd = 0.0
risk_reduction_h  = 0.1
risk_reduction_d  = 0.05

# Direct death probability
θᵍ = Float64.(reduce(hcat, [θᵍ, θᵍ * risk_reduction_dd, θᵍ * risk_reduction_dd]) )
# Hospitalization probability
γᵍ = Float64.(reduce(hcat, [γᵍ, γᵍ * risk_reduction_h, γᵍ * risk_reduction_h]) )
# Fatality probability in ICU
ωᵍ = Float64.(reduce(hcat, [ωᵍ, ωᵍ * risk_reduction_d, ωᵍ * risk_reduction_d]) )

epi_params = Epidemic_Params(βᴵ, βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ, 
                             ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ, 
                             Λ, Γ, rᵥ, kᵥ, 
                             G, M, T)

## Test population parameters
@testset "Epidemic_Params" begin

    # @test epi_params.βᴵ[1] == βᴵ
    # @test epi_params.βᴬ[1] == βᴬ
    # @test epi_params.ηᵍ == ηᵍ
    # @test epi_params.αᵍ == αᵍ
    # @test epi_params.μᵍ == μᵍ
    # @test epi_params.θᵍ == θᵍ
    # @test epi_params.γᵍ == γᵍ
    # @test epi_params.ζᵍ == ζᵍ
    # @test epi_params.λᵍ == λᵍ
    # @test epi_params.ωᵍ == ωᵍ
    # @test epi_params.ψᵍ == ψᵍ
    # @test epi_params.χᵍ == χᵍ
    # @test epi_params.T == T

    @test size(epi_params.ρˢᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴱᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴬᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴵᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴾᴴᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴾᴰᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴴᴿᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴴᴰᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴰᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴿᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.CHᵢᵍᵥ) == (G, M, T, V)
    @test size(epi_params.Qᵢᵍ)   == (G, M, T)

end


## -----------------------------------------------------------------------------
## Loading required data
## -----------------------------------------------------------------------------

data_path       = "data"
config_fname    = joinpath(data_path, "config_test.json")
config          = JSON.parsefile(config_fname)

data_dict       = config["data"]
simulation_dict = config["simulation"]
pop_params_dict = config["population_params"]

# Reading simulation start and end dates
first_day = Date(simulation_dict["first_day_simulation"])
last_day  = Date(simulation_dict["last_day_simulation"])

# Converting dates to time steps
T = (last_day - first_day).value + 1

# Array with time coordinates (dates)
T_coords  = string.(collect(first_day:last_day))

# Loading metapopulation patches info (surface, label, population by age)
metapop_data_filename = joinpath(data_path, data_dict["metapopulation_data_filename"])
metapop_df = CSV.read(metapop_data_filename, DataFrame, types=Dict(:id => String))

# Loading mobility network
mobility_matrix_filename = joinpath(data_path, data_dict["mobility_matrix_filename"])
network_df  = CSV.read(mobility_matrix_filename, DataFrame)

# Metapopulations patches coordinates (labels)
M_coords = map(String,metapop_df[:, "id"])
M = length(M_coords)

# Coordinates for each age strata (labels)
G_coords = map(String, pop_params_dict["age_labels"])
G = length(G_coords)


## -----------------------------------------------------------------------------
## init_pop_param_struct
## -----------------------------------------------------------------------------


population = init_pop_param_struct(G, M, G_coords, pop_params_dict, metapop_df, network_df)
@testset "init_pop_param_struct" begin
    ## POPULATION PARAMETERS
    @test population.G            == G
    @test size(population.C)      == (G, G)
    @test length(population.kᵍ)   == G
    @test length(population.kᵍ_h) == G
    @test length(population.kᵍ_w) == G
    @test length(population.Nᵍ)   == G
    @test length(population.zᵍ)   == G
end



epi_params_dict = config["epidemic_params"]
epi_params      = init_epi_parameters_struct(G, M, T, G_coords, epi_params_dict)
@testset "init_epi_parameters_struct" begin
    # @test epi_params.βᴵ[1] == βᴵ
    # @test epi_params.βᴬ[1] == βᴬ
    @test epi_params.ηᵍ == epi_params_dict["ηᵍ"]
    # @test epi_params.αᵍ == αᵍ
    # @test epi_params.μᵍ == μᵍ
    # @test epi_params.θᵍ == θᵍ
    # @test epi_params.γᵍ == γᵍ
    # @test epi_params.ζᵍ == ζᵍ
    # @test epi_params.λᵍ == λᵍ
    # @test epi_params.ωᵍ == ωᵍ
    # @test epi_params.ψᵍ == ψᵍ
    # @test epi_params.χᵍ == χᵍ
    # @test epi_params.T == T

    @test size(epi_params.ρˢᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴱᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴬᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴵᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴾᴴᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴾᴰᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴴᴿᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴴᴰᵍᵥ) == (G, M, T, V)
    @test size(epi_params.ρᴰᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.ρᴿᵍᵥ)  == (G, M, T, V)
    @test size(epi_params.CHᵢᵍᵥ) == (G, M, T, V)
    @test size(epi_params.Qᵢᵍ)   == (G, M, T)

end



#########################################################
# Containment measures
#########################################################

npi_params_dict = config["NPI"]
# Daily Mobility reduction
kappa0_filename = get(data_dict, "kappa0_filename", nothing)
npi_params = init_NPI_parameters_struct(data_path, npi_params_dict, kappa0_filename, first_day)

#########################################################
# Vaccination parameters
#########################################################

vac_params_dict = config["vaccination"]
start_vacc = vac_params_dict["start_vacc"]
dur_vacc   = vac_params_dict["dur_vacc"]
end_vacc   = start_vacc + dur_vacc

# total vaccinations per age strata
total_population = sum(population.nᵢᵍ)
ϵᵍ = vac_params_dict["ϵᵍ"] * round( total_population * vac_params_dict["percentage_of_vacc_per_day"] )
tᵛs = [start_vacc, end_vacc, T]
ϵᵍs = ϵᵍ .* [0  Int(vac_params_dict["are_there_vaccines"])  0] 


Sᵛ₀ = zeros(Float64, G, M)
E₀  = zeros(Float64, G, M)
A₀  = zeros(Float64, G, M)
I₀  = zeros(Float64, G, M)
H₀  = zeros(Float64, G, M)
R₀  = zeros(Float64, G, M)

A₀[2, 1, 1] = 2.0
A₀[3, 1, 1] = 1.0
I₀[2, 1, 1] = 1.0

set_initial_conditions!(epi_params, population, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)


## Test Initialization of the epidemics
@testset "set_initial_conditions" begin

    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍᵥ[:, :, 1, :]) ≈
        population.N - (sum(E₀) + sum(A₀) + sum(I₀))  atol=0.0001
    
    @test sum(population.nᵢᵍ .* epi_params.ρˢᵍᵥ[:, :, 2, :]) ≈ sum(Sᵛ₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴱᵍᵥ[:, :, 1, :]) ≈ sum(E₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴬᵍᵥ[:, :, 1, :]) ≈ sum(A₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴵᵍᵥ[:, :, 1, :]) ≈ sum(I₀)  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴴᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴾᴰᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴿᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴴᴰᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴰᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.ρᴿᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001
    @test sum(population.nᵢᵍ .* epi_params.CHᵢᵍᵥ[:, :, 1, :]) ≈ 0.0  atol=0.0001

end


# @test size(initial_compartments) == (G, M, epi_params.V, epi_params.NumComps)
# set_compartments!(epi_params, population, initial_compartments)
# run_epidemic_spreading_mmca!(epi_params, population, npi_params, tᵛs, ϵᵍs; verbose = true )