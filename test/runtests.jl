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



@testset "Parse_config" begin
    data_path = "data"
    config_fname = joinpath(data_path, "config_test.json")
    config = JSON.parsefile(config_fname)

    simulation_dict = config["simulation"]
    data_dict       = config["data"]
    epi_params_dict = config["epidemic_params"]
    pop_params_dict = config["population_params"]
    vac_params_dict = config["vaccination"]
    npi_params_dict = config["NPI"]

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

    ## POPULATION PARAMETERS
    population       = init_pop_param_struct(G, M, G_coords, pop_params_dict, metapop_df, network_df)
    
    ## EPIDEMIC PARAMETERS 
    epi_params       = init_epi_parameters_struct(G, M, T, G_coords, epi_params_dict)

    # @test size(initial_compartments) == (G, M, epi_params.V, epi_params.NumComps)
    @test population.G            == G
    @test size(population.C)      == (G, G)
    @test length(population.kᵍ)   == G
    @test length(population.kᵍ_h) == G
    @test length(population.kᵍ_w) == G
    @test length(population.Nᵍ)   == G
    @test length(population.zᵍ)   == G
    # @test size(population.nᵢᵍ_eff) == (G, M)
    # @test size(population.normᵍ) == (G, M)
    # @test size(population.mobilityᵍ) == (G, length(Rᵢⱼ))

end
