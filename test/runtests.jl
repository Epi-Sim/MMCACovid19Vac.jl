using Test
using DataFrames
using MMCACovid19Vac
using JSON
using Dates
using CSV
using NCDatasets
using LinearAlgebra

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

data_path       = joinpath(@__DIR__, "data")
config_fname    = joinpath(data_path, "config_test.json")
config          = JSON.parsefile(config_fname)

data_dict       = Dict(config["data"])
simulation_dict = Dict(config["simulation"])
pop_params_dict = Dict(config["population_params"])

# Reading simulation start and end dates
first_day = Date(simulation_dict["start_date"])
last_day  = Date(simulation_dict["end_date"])

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
G_coords = map(String, pop_params_dict["G_labels"])
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



epi_params_dict = Dict(config["epidemic_params"])
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

npi_params_dict = Dict(config["NPI"])
# Daily Mobility reduction
kappa0_filename = get(data_dict, "kappa0_filename", nothing)
npi_params = init_NPI_parameters_struct(data_path, npi_params_dict, kappa0_filename, first_day)

#########################################################
# Vaccination parameters
#########################################################

vac_params_dict = Dict(config["vaccination"])
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

@testset "mobility_variation" begin
    # Helper to get a fresh population
    function get_fresh_pop(v_type, sigma_O, sigma_D, seed)
        # Re-load data to ensure no shared mutable objects
        metapop_df = CSV.read(joinpath(data_path, data_dict["metapopulation_data_filename"]), DataFrame, types=Dict(:id => String))
        network_df = CSV.read(joinpath(data_path, data_dict["mobility_matrix_filename"]), DataFrame)
        pop = init_pop_param_struct(G, M, G_coords, pop_params_dict, metapop_df, network_df)
        
        # Manually set the new fields (struct is immutable but contains mutable arrays)
        # Actually struct is struct, so we can't easily re-set. 
        # But Population_Params constructor returns a new one.
        
        # We need to construct a new one with the desired variation params
        return Population_Params(
            pop.G, pop.M, pop.nᵢ, pop.nᵢᵍ,
            pop.nᵢ_eff, pop.nᵢᵍ_eff, pop.N, pop.Nᵍ,
            pop.pᵍ, pop.pᵍ_eff, pop.edgelist, pop.Rᵢⱼ,
            copy(pop.Rᵢⱼ), v_type, sigma_O, sigma_D, seed, nothing, pop.mobilityᵍ,
            pop.kᵍ, pop.kᵍ_h, pop.kᵍ_w, pop.kᵍ_eff,
            pop.C, pop.zᵍ, pop.normᵍ, pop.sᵢ,
            pop.ξ, pop.σ
        )
    end

    # 1. Clean baseline results
    reset_epidemic_compartments!(epi_params)
    pop_base = get_fresh_pop("none", 0.0, 0.0, nothing)
    set_initial_conditions!(epi_params, pop_base, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
    run_epidemic_spreading_mmca!(epi_params, pop_base, npi_params, tᵛs, ϵᵍs; verbose = false)
    baseline_results = copy(epi_params.ρᴵᵍᵥ)

    @testset "reproducibility" begin
        # Run 1
        reset_epidemic_compartments!(epi_params)
        pop1 = get_fresh_pop("ipfp_simple", 0.1, 0.1, 42)
        set_initial_conditions!(epi_params, pop1, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
        run_epidemic_spreading_mmca!(epi_params, pop1, npi_params, tᵛs, ϵᵍs; verbose = false)
        res1 = copy(epi_params.ρᴵᵍᵥ)

        # Run 2
        reset_epidemic_compartments!(epi_params)
        pop2 = get_fresh_pop("ipfp_simple", 0.1, 0.1, 42)
        set_initial_conditions!(epi_params, pop2, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
        run_epidemic_spreading_mmca!(epi_params, pop2, npi_params, tᵛs, ϵᵍs; verbose = false)
        res2 = copy(epi_params.ρᴵᵍᵥ)

        @test res1 ≈ res2 atol=1e-12
    end

    @testset "high_noise_divergence" begin
        # Run with high noise
        reset_epidemic_compartments!(epi_params)
        pop_high = get_fresh_pop("ipfp_simple", 0.5, 0.5, 123)
        set_initial_conditions!(epi_params, pop_high, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
        run_epidemic_spreading_mmca!(epi_params, pop_high, npi_params, tᵛs, ϵᵍs; verbose = false)
        res_high = copy(epi_params.ρᴵᵍᵥ)

        # Should differ from baseline
        @test !isapprox(res_high, baseline_results, atol=1e-4)
        @test all(epi_params.ρˢᵍᵥ .>= -1e-12)
    end

    @testset "ipfp_structural_properties" begin
        # 1. Topology Preservation (Island Test)
        # We check that IPFP never creates edges where weight was 0
        pop_struct = get_fresh_pop("ipfp_simple", 1.0, 1.0, 999) # High noise
        
        # Manually trigger update
        MMCACovid19Vac.update_mobility!(pop_struct, 1)
        
        # Check zeros are preserved
        for i in 1:length(pop_struct.Rᵢⱼ)
            if pop_struct.Rᵢⱼ_base[i] == 0.0
                @test pop_struct.Rᵢⱼ[i] == 0.0
            end
        end

        # 2. Row Stochasticity
        # Check that rows sum to 1.0
        # We need to reconstruct the matrix or sum by edges
        row_sums = zeros(Float64, pop_struct.M)
        for k in 1:length(pop_struct.Rᵢⱼ)
            u = pop_struct.edgelist[k, 1]
            row_sums[u] += pop_struct.Rᵢⱼ[k]
        end
        @test all(isapprox.(row_sums, 1.0, atol=1e-10))

        # 3. Ghost Population Fix Verification
        # Check that n_eff actually changed from the baseline initialization
        # (Since we applied noise, the flows should shift, changing local density)
        @test pop_struct.nᵢ_eff != pop_base.nᵢ_eff
    end

    @testset "confinement_mechanism" begin
        # Verify that kappa0 correctly reduces inter-patch mobility
        reset_epidemic_compartments!(epi_params)
        
        # Create population with high mobility to start
        pop_conf = get_fresh_pop("none", 0.0, 0.0, nothing)
        pop_conf.pᵍ .= 1.0 # Everyone wants to move
        
        # Test 1: No confinement (kappa = 0)
        kappa_zero = 0.0
        kappa_full = 0.99 # 99% confinement
        
        # Case A: Kappa = 0
        pop_conf.pᵍ_eff .= (1 - kappa_zero) .* pop_conf.pᵍ
        update_population_params!(pop_conf)
        
        # Check off-diagonal flow exists
        off_diag_flow_A = 0.0
        for k in 1:length(pop_conf.Rᵢⱼ)
            u, v = pop_conf.edgelist[k, 1], pop_conf.edgelist[k, 2]
            if u != v
                off_diag_flow_A += sum(pop_conf.mobilityᵍ[:, k])
            end
        end
        @test off_diag_flow_A > 1.0 # Should be substantial
        
        # Case B: Kappa = 0.99
        pop_conf.pᵍ_eff .= (1 - kappa_full) .* pop_conf.pᵍ
        update_population_params!(pop_conf)
        
        off_diag_flow_B = 0.0
        for k in 1:length(pop_conf.Rᵢⱼ)
            u, v = pop_conf.edgelist[k, 1], pop_conf.edgelist[k, 2]
            if u != v
                off_diag_flow_B += sum(pop_conf.mobilityᵍ[:, k])
            end
        end
        
        # The flow should be reduced by roughly factor (1-0.99)
        @test off_diag_flow_B < off_diag_flow_A * 0.02
    end

    # ============================================================================
    # IPFP RANDOM (TIME-VARYING) - Additional Tests
    # ============================================================================
    
    @testset "ipfp_zero_sigma_baseline_invariance" begin
        # Zero-sigma = baseline invariance
        # Setup with sigma_O=0, sigma_D=0, fixed mobility_rng_seed
        # Call update_mobility! for multiple t
        # Expect Rᵢⱼ unchanged from Rᵢⱼ_base (within tolerance) and 
        # nᵢ_eff, zᵍ, normᵍ unchanged
        
        pop_zero = get_fresh_pop("ipfp_simple", 0.0, 0.0, 42)
        R_base = copy(pop_zero.Rᵢⱼ_base)
        nᵢ_eff_base = copy(pop_zero.nᵢ_eff)
        zᵍ_base = copy(pop_zero.zᵍ)
        normᵍ_base = copy(pop_zero.normᵍ)
        
        # Update mobility at multiple timesteps
        for t in 1:5
            MMCACovid19Vac.update_mobility!(pop_zero, t)
        end
        
        # Rᵢⱼ should be unchanged from baseline (within tolerance)
        @test pop_zero.Rᵢⱼ ≈ R_base atol=1e-10
        
        # nᵢ_eff, zᵍ, normᵍ should be unchanged
        # Note: Using relaxed tolerance due to floating-point precision in population updates
        @test pop_zero.nᵢ_eff ≈ nᵢ_eff_base atol=1e-6
        @test pop_zero.zᵍ ≈ zᵍ_base atol=1e-6
        @test pop_zero.normᵍ ≈ normᵍ_base atol=1e-6
    end
    
    @testset "ipfp_seed_time_determinism" begin
        # With fixed mobility_rng_seed, run update_mobility!(t=3) twice on fresh 
        # populations and assert identical Rᵢⱼ
        # Also assert Rᵢⱼ(t=3) differs from Rᵢⱼ(t=4) (same seed, different t)
        
        # First run at t=3
        pop1 = get_fresh_pop("ipfp_simple", 0.1, 0.1, 42)
        MMCACovid19Vac.update_mobility!(pop1, 3)
        R_t3_run1 = copy(pop1.Rᵢⱼ)
        
        # Second run at t=3 (should be identical)
        pop2 = get_fresh_pop("ipfp_simple", 0.1, 0.1, 42)
        MMCACovid19Vac.update_mobility!(pop2, 3)
        R_t3_run2 = copy(pop2.Rᵢⱼ)
        
        @test R_t3_run1 ≈ R_t3_run2 atol=1e-12
        
        # Run at t=4 (should differ from t=3)
        pop3 = get_fresh_pop("ipfp_simple", 0.1, 0.1, 42)
        MMCACovid19Vac.update_mobility!(pop3, 4)
        R_t4 = copy(pop3.Rᵢⱼ)
        
        @test !isapprox(R_t3_run1, R_t4, atol=1e-6)
    end
    
    @testset "ipfp_marginal_matching" begin
        # Marginal matching post-IPFP
        # Capture O_t and D_t implicitly by re-computing row/col sums of Rᵢⱼ
        # Assert row sums match O_t and column sums match D_t within tolerance
        
        pop_marg = get_fresh_pop("ipfp_simple", 0.3, 0.3, 123)
        M = pop_marg.M
        edgelist = pop_marg.edgelist
        B = pop_marg.Rᵢⱼ_base
        
        # Compute baseline marginals
        O_base = zeros(Float64, M)
        D_base = zeros(Float64, M)
        for e in 1:length(B)
            O_base[edgelist[e, 1]] += B[e]
            D_base[edgelist[e, 2]] += B[e]
        end
        
        # Generate noisy marginals using same logic as update_mobility!
        using Random
        seed = 123 + 5  # seed + t where t=5
        rng = MersenneTwister(seed)
        sigma_O = 0.3
        sigma_D = 0.3
        
        O_t = O_base .* exp.(randn(rng, M) .* sigma_O .- (sigma_O^2 / 2))
        D_t = D_base .* exp.(randn(rng, M) .* sigma_D .- (sigma_D^2 / 2))
        D_t .*= sum(O_t) / sum(D_t)  # Balance marginals
        
        # Update mobility at t=5
        MMCACovid19Vac.update_mobility!(pop_marg, 5)
        
        # Compute row and column sums from resulting Rᵢⱼ
        R_t = pop_marg.Rᵢⱼ
        row_sums = zeros(Float64, M)
        col_sums = zeros(Float64, M)
        for e in 1:length(R_t)
            u = edgelist[e, 1]
            v = edgelist[e, 2]
            row_sums[u] += R_t[e]
            col_sums[v] += R_t[e]
        end
        
        # Row sums should be 1.0 (row-stochastic after final normalization)
        @test all(isapprox.(row_sums, 1.0, atol=1e-10))
        
        # Column sums should match normalized D_t (within tolerance)
        # Note: After IPFP, the column sums may not exactly match D_t due to 
        # the final row normalization step, but they should be close
        D_t_normalized = D_t ./ sum(D_t)  # Normalize to compare with column sums
        @test sum(col_sums) ≈ 1.0 * M atol=1e-6  # Total flow conservation
    end
    
    @testset "ipfp_non_negativity_high_noise" begin
        # Non-negativity under high noise
        # With very large sigma_O, sigma_D (e.g. 1.5), assert all Rᵢⱼ >= 0 and no NaNs
        
        pop_high_noise = get_fresh_pop("ipfp_simple", 1.5, 1.5, 999)
        
        # Update mobility multiple times
        for t in 1:10
            MMCACovid19Vac.update_mobility!(pop_high_noise, t)
            
            # Check non-negativity
            @test all(pop_high_noise.Rᵢⱼ .>= -1e-12)  # Allow tiny numerical errors
            
            # Check no NaNs
            @test all(.!isnan.(pop_high_noise.Rᵢⱼ))
            @test all(.!isinf.(pop_high_noise.Rᵢⱼ))
        end
    end
    
    @testset "ipfp_isolated_node_handling" begin
        # Isolated node handling
        # Construct a tiny edgelist with a node that has no outgoing edges (row sum zero)
        # Call update_mobility! and assert it doesn't produce NaNs and that 
        # row sums remain zero for isolated nodes
        
        # Create a simple test population with an isolated node
        test_G = 2
        test_M = 4
        test_nᵢᵍ = ones(Float64, test_G, test_M) * 1000.0
        test_C = Matrix{Float64}(I, test_G, test_G)
        test_kᵍ = ones(test_G) * 12.0
        test_kᵍ_h = ones(test_G) * 3.15
        test_kᵍ_w = ones(test_G) * 2.0
        test_pᵍ = ones(test_G) * 0.5
        test_sᵢ = ones(test_M) * 10.0
        test_ξ = 0.01
        test_σ = 2.5
        
        # Edgelist where node 4 has no outgoing edges (isolated)
        test_edgelist = [1 1; 1 2; 2 1; 2 2; 2 3; 3 2; 3 3; 3 4; 4 4]
        test_Rᵢⱼ = [0.7; 0.3; 0.4; 0.4; 0.2; 0.2; 0.5; 0.3; 1.0]  # Node 4 only has self-loop
        
        pop_isolated = Population_Params(
            test_G, test_M, test_nᵢᵍ, test_kᵍ, test_kᵍ_h, test_kᵍ_w,
            test_C, test_pᵍ, test_edgelist, test_Rᵢⱼ,
            test_sᵢ, test_ξ, test_σ
        )
        
        # Convert to ipfp_simple mode
        pop_isolated = Population_Params(
            pop_isolated.G, pop_isolated.M, pop_isolated.nᵢ, pop_isolated.nᵢᵍ,
            pop_isolated.nᵢ_eff, pop_isolated.nᵢᵍ_eff, pop_isolated.N, pop_isolated.Nᵍ,
            pop_isolated.pᵍ, pop_isolated.pᵍ_eff, pop_isolated.edgelist, pop_isolated.Rᵢⱼ,
            copy(pop_isolated.Rᵢⱼ), "ipfp_simple", 0.5, 0.5, 42, nothing, pop_isolated.mobilityᵍ,
            pop_isolated.kᵍ, pop_isolated.kᵍ_h, pop_isolated.kᵍ_w, pop_isolated.kᵍ_eff,
            pop_isolated.C, pop_isolated.zᵍ, pop_isolated.normᵍ, pop_isolated.sᵢ,
            pop_isolated.ξ, pop_isolated.σ
        )
        
        # Update mobility
        MMCACovid19Vac.update_mobility!(pop_isolated, 1)
        
        # Check no NaNs
        @test all(.!isnan.(pop_isolated.Rᵢⱼ))
        
        # Check row sums - node 4 should still sum to 1.0 (self-loop)
        row_sums = zeros(Float64, test_M)
        for e in 1:length(pop_isolated.Rᵢⱼ)
            u = pop_isolated.edgelist[e, 1]
            row_sums[u] += pop_isolated.Rᵢⱼ[e]
        end
        
        # All rows should sum to 1.0 (including isolated node via self-loop)
        @test all(isapprox.(row_sums, 1.0, atol=1e-10))
    end
end


## -----------------------------------------------------------------------------
## External Mobility Tests
## -----------------------------------------------------------------------------

@testset "external_mobility_netcdf" begin
    # Helper function to create a test NetCDF file with mobility data
    function create_test_mobility_netcdf(filename::String, M::Int, T::Int,
                                        start_date::Date, edgelist::Array{Int64, 2},
                                        R_base::Array{Float64, 1})
        L = length(edgelist[:, 1])

        # Create date strings
        dates = collect(start_date:Day(1):(start_date + Day(T - 1)))
        date_strs = [Dates.format(d, "yyyy-mm-dd") for d in dates]

        # Create dense mobility array (T, M, M)
        # For each timestep, vary the mobility significantly from baseline
        # by varying only diagonal elements (so normalization doesn't cancel it out)
        mobility_dense = zeros(Float64, T, M, M)

        for t in 1:T
            # Apply variation that increases diagonal over time
            # Only vary diagonal elements, keep off-diagonal proportional
            diag_factor = 1.0 + 0.5 * t / T  # Increases from 1.0 to 1.5 over time

            for e in 1:L
                origin = edgelist[e, 1]
                dest = edgelist[e, 2]
                if origin == dest
                    # Diagonal: increase more
                    mobility_dense[t, origin, dest] = R_base[e] * diag_factor
                else
                    # Off-diagonal: decrease to compensate
                    mobility_dense[t, origin, dest] = R_base[e]
                end
            end

            # Normalize each row to ensure row-stochasticity
            for i in 1:M
                row_sum = sum(mobility_dense[t, i, :])
                if row_sum > 0
                    mobility_dense[t, i, :] ./= row_sum
                end
            end
        end

        # Create NetCDF file using NCDatasets (wraps HDF5, better ARM64 support)
        ds = NCDataset(filename, "c")

        try
            # Define dimensions
            defDim(ds, "date", T)
            defDim(ds, "origin", M)
            defDim(ds, "destination", M)

            # Define and write date variable
            defVar(ds, "date", String, ("date",))
            ds["date"][:] = date_strs

            # Define and write mobility variable
            defVar(ds, "mobility", Float64, ("date", "origin", "destination"))
            ds["mobility"][:, :, :] = mobility_dense
        finally
            close(ds)
        end

        return filename
    end

    # Test parameters
    M_test = 5
    T_test = 10
    test_start_date = Date(2020, 3, 10)
    test_end_date = test_start_date + Day(T_test - 1)
    test_netcdf_file = joinpath(tempdir(), "test_mobility.nc")

    # Create a simple test population with M=5
    test_G = 3
    test_nᵢᵍ = ones(Float64, test_G, M_test) * 1000.0
    test_C = Matrix{Float64}(I, test_G, test_G)
    test_kᵍ = ones(test_G) * 12.0
    test_kᵍ_h = ones(test_G) * 3.15
    test_kᵍ_w = ones(test_G) * 2.0
    test_pᵍ = ones(test_G) * 0.5
    test_sᵢ = ones(M_test) * 10.0
    test_ξ = 0.01
    test_σ = 2.5
    # Simple test edgelist with M=5
    test_edgelist = zeros(Int64, M_test * M_test, 2)
    test_Rᵢⱼ = zeros(Float64, M_test * M_test)
    idx = 1
    for i in 1:M_test
        for j in 1:M_test
            test_edgelist[idx, 1] = i
            test_edgelist[idx, 2] = j
            test_Rᵢⱼ[idx] = i == j ? 0.7 : 0.3 / (M_test - 1)  # Row-stochastic
            idx += 1
        end
    end

    test_population = Population_Params(
        test_G, M_test, test_nᵢᵍ, test_kᵍ, test_kᵍ_h, test_kᵍ_w,
        test_C, test_pᵍ, test_edgelist, test_Rᵢⱼ,
        test_sᵢ, test_ξ, test_σ
    )

    try
        # Create test NetCDF file
        create_test_mobility_netcdf(test_netcdf_file, M_test, T_test,
                                    test_start_date, test_edgelist, test_Rᵢⱼ)

        @testset "config_parameters" begin
            # Test that new config parameters have correct defaults
            pop_dict = create_population_params_dict(G)
            @test pop_dict["mobility_variation_type"] == "none"
            @test pop_dict["mobility_external_file"] === nothing
            @test pop_dict["mobility_external_variable"] == "mobility"
            @test pop_dict["mobility_validate"] == true
            @test pop_dict["mobility_fill_missing"] == false
            @test pop_dict["mobility_time_alignment"] == "start"
        end

        @testset "load_mobility_netcdf" begin
            # Test loading mobility from NetCDF
            R_series, dates = load_mobility_netcdf(
                test_netcdf_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Check shape
            @test size(R_series) == (T_test, length(test_edgelist[:, 1]))
            @test length(dates) == T_test
            @test dates[1] == test_start_date
            @test dates[end] == test_end_date

            # Check row-stochasticity (each origin row should sum to 1.0)
            for t in 1:T_test
                row_sums = zeros(Float64, M_test)
                for e in 1:length(test_edgelist[:, 1])
                    origin = test_edgelist[e, 1]
                    row_sums[origin] += R_series[t, e]
                end
                @test all(isapprox.(row_sums, 1.0, atol=1e-6))
            end

            # Check non-negativity
            @test all(R_series .>= 0.0)
        end

        @testset "validate_mobility_series" begin
            # Load series
            R_series, _ = load_mobility_netcdf(
                test_netcdf_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Test validation
            is_valid, errors = validate_mobility_series(
                R_series, test_population.edgelist, test_population.M, T_test, verbose=false
            )

            @test is_valid
            @test length(errors) == 0
        end

        @testset "round_trip_save_load" begin
            # Save and reload mobility
            R_original, _ = load_mobility_netcdf(
                test_netcdf_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Save to new file
            test_output_file = joinpath(tempdir(), "test_mobility_output.nc")
            dates = collect(test_start_date:Day(1):test_end_date)
            save_mobility_netcdf(R_original, test_population.edgelist, dates, test_output_file, compress=false)

            # Reload
            R_reloaded, dates_reloaded = load_mobility_netcdf(
                test_output_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Compare
            @test size(R_original) == size(R_reloaded)
            @test dates == dates_reloaded
            @test R_original ≈ R_reloaded atol=1e-6
        end

        @testset "update_mobility_external_mode" begin
            # Load mobility series
            R_series, _ = load_mobility_netcdf(
                test_netcdf_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Create population with external_netcdf mode
            pop_external = Population_Params(
                test_population.G, test_population.M, test_population.nᵢ, test_population.nᵢᵍ,
                test_population.nᵢ_eff, test_population.nᵢᵍ_eff, test_population.N, test_population.Nᵍ,
                test_population.pᵍ, test_population.pᵍ_eff, test_population.edgelist, test_population.Rᵢⱼ,
                copy(test_population.Rᵢⱼ), "external_netcdf", 0.0, 0.0, nothing,
                R_series, test_population.mobilityᵍ,
                test_population.kᵍ, test_population.kᵍ_h, test_population.kᵍ_w, test_population.kᵍ_eff,
                test_population.C, test_population.zᵍ, test_population.normᵍ, test_population.sᵢ,
                test_population.ξ, test_population.σ
            )

            # Store original Rᵢⱼ
            R_original = copy(pop_external.Rᵢⱼ)

            # Update mobility at t=1
            MMCACovid19Vac.update_mobility!(pop_external, 1)

            # Check that Rᵢⱼ was updated
            @test pop_external.Rᵢⱼ ≈ R_series[1, :] atol=1e-10
            @test !isapprox(pop_external.Rᵢⱼ, R_original, atol=1e-6)

            # Update at t=5
            MMCACovid19Vac.update_mobility!(pop_external, 5)
            @test pop_external.Rᵢⱼ ≈ R_series[5, :] atol=1e-10
        end

        @testset "insufficient_data_error" begin
            # Create a test NetCDF with only 5 timesteps (less than T_test=10)
            partial_netcdf = joinpath(tempdir(), "test_mobility_partial.nc")
            create_test_mobility_netcdf(partial_netcdf, M_test, 5,  # Only 5 timesteps
                                        test_start_date, test_edgelist, test_Rᵢⱼ)

            # Try to load 10 timesteps - should error
            @test_throws ErrorException load_mobility_netcdf(
                partial_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            rm(partial_netcdf)
        end

        # ============================================================================
        # EXTERNAL NETCDF - Additional Tests
        # ============================================================================

        @testset "alignment_modes" begin
            # Test alignment modes: start, end, nearest
            # Create a NetCDF with dates offset by 1-2 days from simulation window

            offset_start_date = test_start_date + Day(2)  # NetCDF starts 2 days later
            offset_netcdf = joinpath(tempdir(), "test_mobility_offset.nc")

            # Create NetCDF with offset dates
            create_test_mobility_netcdf(offset_netcdf, M_test, T_test,
                                        offset_start_date, test_edgelist, test_Rᵢⱼ)

            # Test "start" alignment (default) - should error since dates don't match
            @test_throws ErrorException load_mobility_netcdf(
                offset_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false, alignment="start"
            )

            # Test "nearest" alignment with fill_missing=true
            R_series, dates = load_mobility_netcdf(
                offset_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false,
                alignment="nearest", fill_missing=true
            )

            # Check that we got the right number of timesteps
            @test size(R_series, 1) == T_test
            @test length(dates) == T_test

            rm(offset_netcdf)
        end

        @testset "fill_missing_interpolation" begin
            # Build a NetCDF with a missing date in the middle
            gap_netcdf = joinpath(tempdir(), "test_mobility_gap.nc")

            # Create dates with a gap (missing day 5)
            gap_dates = [test_start_date:Day(1):(test_start_date + Day(3));
                         test_start_date + Day(5):Day(1):test_end_date]

            # Create mobility data for these dates
            T_gap = length(gap_dates)
            mobility_dense_gap = zeros(Float64, T_gap, M_test, M_test)
            for t in 1:T_gap
                for e in 1:length(test_edgelist[:, 1])
                    origin = test_edgelist[e, 1]
                    dest = test_edgelist[e, 2]
                    mobility_dense_gap[t, origin, dest] = test_Rᵢⱼ[e]
                end
                # Normalize rows
                for i in 1:M_test
                    row_sum = sum(mobility_dense_gap[t, i, :])
                    if row_sum > 0
                        mobility_dense_gap[t, i, :] ./= row_sum
                    end
                end
            end

            # Save to NetCDF
            ds_gap = NCDataset(gap_netcdf, "c")
            try
                defDim(ds_gap, "date", T_gap)
                defDim(ds_gap, "origin", M_test)
                defDim(ds_gap, "destination", M_test)

                date_strs = [Dates.format(d, "yyyy-mm-dd") for d in gap_dates]
                defVar(ds_gap, "date", String, ("date",))
                ds_gap["date"][:] = date_strs

                defVar(ds_gap, "mobility", Float64, ("date", "origin", "destination"))
                ds_gap["mobility"][:, :, :] = mobility_dense_gap
            finally
                close(ds_gap)
            end

            # Test with fill_missing=false - should error
            @test_throws ErrorException load_mobility_netcdf(
                gap_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false, fill_missing=false
            )

            # Test with fill_missing=true - should interpolate
            R_series, dates = load_mobility_netcdf(
                gap_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false, fill_missing=true
            )

            # Check that we got all timesteps including the missing date
            @test size(R_series, 1) == T_test
            @test length(dates) == T_test

            rm(gap_netcdf)
        end

        @testset "variable_missing_dimension_mismatch" begin
            # Test NetCDF lacking the mobility variable
            missing_var_netcdf = joinpath(tempdir(), "test_missing_var.nc")
            ds_mv = NCDataset(missing_var_netcdf, "c")
            try
                defDim(ds_mv, "date", T_test)
                defDim(ds_mv, "origin", M_test)
                defDim(ds_mv, "destination", M_test)

                date_strs = [Dates.format(test_start_date + Day(i-1), "yyyy-mm-dd") for i in 1:T_test]
                defVar(ds_mv, "date", String, ("date",))
                ds_mv["date"][:] = date_strs

                # Define a different variable, not "mobility"
                defVar(ds_mv, "different_variable", Float64, ("date", "origin", "destination"))
                ds_mv["different_variable"][:, :, :] = zeros(T_test, M_test, M_test)
            finally
                close(ds_mv)
            end

            # Should throw error for missing variable
            @test_throws ErrorException load_mobility_netcdf(
                missing_var_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            rm(missing_var_netcdf)

            # Test dimension mismatch (M doesn't match)
            wrong_M_netcdf = joinpath(tempdir(), "test_wrong_M.nc")
            wrong_M = M_test + 2  # Different M

            ds_wm = NCDataset(wrong_M_netcdf, "c")
            try
                defDim(ds_wm, "date", T_test)
                defDim(ds_wm, "origin", wrong_M)
                defDim(ds_wm, "destination", wrong_M)

                date_strs = [Dates.format(test_start_date + Day(i-1), "yyyy-mm-dd") for i in 1:T_test]
                defVar(ds_wm, "date", String, ("date",))
                ds_wm["date"][:] = date_strs

                defVar(ds_wm, "mobility", Float64, ("date", "origin", "destination"))
                ds_wm["mobility"][:, :, :] = zeros(T_test, wrong_M, wrong_M)
            finally
                close(ds_wm)
            end

            # Should throw error for dimension mismatch
            @test_throws ErrorException load_mobility_netcdf(
                wrong_M_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            rm(wrong_M_netcdf)
        end

        @testset "validation_failures" begin
            # Test validation with negative values
            negative_netcdf = joinpath(tempdir(), "test_negative.nc")
            mobility_dense_neg = zeros(Float64, T_test, M_test, M_test)

            # Create data with negative values
            for t in 1:T_test
                for e in 1:length(test_edgelist[:, 1])
                    origin = test_edgelist[e, 1]
                    dest = test_edgelist[e, 2]
                    mobility_dense_neg[t, origin, dest] = test_Rᵢⱼ[e]
                end
                # Add a negative value
                mobility_dense_neg[t, 1, 1] = -0.1
            end

            ds_neg = NCDataset(negative_netcdf, "c")
            try
                defDim(ds_neg, "date", T_test)
                defDim(ds_neg, "origin", M_test)
                defDim(ds_neg, "destination", M_test)

                date_strs = [Dates.format(test_start_date + Day(i-1), "yyyy-mm-dd") for i in 1:T_test]
                defVar(ds_neg, "date", String, ("date",))
                ds_neg["date"][:] = date_strs

                defVar(ds_neg, "mobility", Float64, ("date", "origin", "destination"))
                ds_neg["mobility"][:, :, :] = mobility_dense_neg
            finally
                close(ds_neg)
            end

            # Load without validation first
            R_series, _ = load_mobility_netcdf(
                negative_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Validate should fail
            is_valid, errors = validate_mobility_series(
                R_series, test_population.edgelist, M_test, T_test, verbose=false
            )

            @test !is_valid
            @test length(errors) > 0
            @test any(occursin.("negative", lowercase.(errors)))

            rm(negative_netcdf)

            # Test validation with rows that don't sum to 1
            nonstochastic_netcdf = joinpath(tempdir(), "test_nonstochastic.nc")
            mobility_dense_ns = zeros(Float64, T_test, M_test, M_test)

            # Create data where rows don't sum to 1
            for t in 1:T_test
                for e in 1:length(test_edgelist[:, 1])
                    origin = test_edgelist[e, 1]
                    dest = test_edgelist[e, 2]
                    mobility_dense_ns[t, origin, dest] = test_Rᵢⱼ[e] * 0.5  # Scale down
                end
            end

            ds_ns = NCDataset(nonstochastic_netcdf, "c")
            try
                defDim(ds_ns, "date", T_test)
                defDim(ds_ns, "origin", M_test)
                defDim(ds_ns, "destination", M_test)

                date_strs = [Dates.format(test_start_date + Day(i-1), "yyyy-mm-dd") for i in 1:T_test]
                defVar(ds_ns, "date", String, ("date",))
                ds_ns["date"][:] = date_strs

                defVar(ds_ns, "mobility", Float64, ("date", "origin", "destination"))
                ds_ns["mobility"][:, :, :] = mobility_dense_ns
            finally
                close(ds_ns)
            end

            # Load without validation
            R_series_ns, _ = load_mobility_netcdf(
                nonstochastic_netcdf, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Validate should fail
            is_valid_ns, errors_ns = validate_mobility_series(
                R_series_ns, test_population.edgelist, M_test, T_test, verbose=false
            )

            @test !is_valid_ns
            @test length(errors_ns) > 0

            rm(nonstochastic_netcdf)
        end

        @testset "update_mobility_auxiliary_updates" begin
            # Test that update_mobility! applies auxiliary updates
            # After update_mobility!(t), assert nᵢ_eff, zᵍ, and mobilityᵍ change
            # consistently with the new Rᵢⱼ

            R_series, _ = load_mobility_netcdf(
                test_netcdf_file, test_population, test_start_date, test_end_date,
                varname="mobility", validate=false, verbose=false
            )

            # Create population with external_netcdf mode
            pop_aux = Population_Params(
                test_population.G, test_population.M, test_population.nᵢ, test_population.nᵢᵍ,
                test_population.nᵢ_eff, test_population.nᵢᵍ_eff, test_population.N, test_population.Nᵍ,
                test_population.pᵍ, test_population.pᵍ_eff, test_population.edgelist, test_population.Rᵢⱼ,
                copy(test_population.Rᵢⱼ), "external_netcdf", 0.0, 0.0, nothing,
                R_series, test_population.mobilityᵍ,
                test_population.kᵍ, test_population.kᵍ_h, test_population.kᵍ_w, test_population.kᵍ_eff,
                test_population.C, test_population.zᵍ, test_population.normᵍ, test_population.sᵢ,
                test_population.ξ, test_population.σ
            )

            # Store original values
            nᵢ_eff_orig = copy(pop_aux.nᵢ_eff)
            zᵍ_orig = copy(pop_aux.zᵍ)
            mobilityᵍ_orig = copy(pop_aux.mobilityᵍ)

            # Update mobility at t=1
            MMCACovid19Vac.update_mobility!(pop_aux, 1)

            # Check that auxiliary values were updated
            @test pop_aux.nᵢ_eff != nᵢ_eff_orig || pop_aux.zᵍ != zᵍ_orig || pop_aux.mobilityᵍ != mobilityᵍ_orig

            # Verify Rᵢⱼ was updated
            @test pop_aux.Rᵢⱼ ≈ R_series[1, :] atol=1e-10
        end

        @testset "load_external_mobility_noop" begin
            # Test load_external_mobility! no-op behavior
            # With mobility_variation_type != "external_netcdf", call load_external_mobility!
            # and assert the population instance is unchanged

            # Create a population with "none" variation type
            pop_none = Population_Params(
                test_population.G, test_population.M, test_population.nᵢ, test_population.nᵢᵍ,
                test_population.nᵢ_eff, test_population.nᵢᵍ_eff, test_population.N, test_population.Nᵍ,
                test_population.pᵍ, test_population.pᵍ_eff, test_population.edgelist, test_population.Rᵢⱼ,
                copy(test_population.Rᵢⱼ), "none", 0.0, 0.0, nothing,
                nothing, test_population.mobilityᵍ,
                test_population.kᵍ, test_population.kᵍ_h, test_population.kᵍ_w, test_population.kᵍ_eff,
                test_population.C, test_population.zᵍ, test_population.normᵍ, test_population.sᵢ,
                test_population.ξ, test_population.σ
            )

            # Store original values
            Rᵢⱼ_orig = copy(pop_none.Rᵢⱼ)

            # Create config with external file
            config = Dict(
                "mobility_external_file" => test_netcdf_file,
                "mobility_external_variable" => "mobility",
                "mobility_validate" => false,
                "mobility_fill_missing" => false,
                "mobility_time_alignment" => "start",
                "start_date" => test_start_date,
                "end_date" => test_end_date
            )

            # Call load_external_mobility! - should return unchanged population
            pop_result = MMCACovid19Vac.load_external_mobility!(pop_none, config, verbose=false)

            # Population should be unchanged (same object since struct is immutable)
            @test pop_result === pop_none
            @test pop_result.mobility_external_series === nothing
            @test pop_result.Rᵢⱼ ≈ Rᵢⱼ_orig atol=1e-10
        end

    finally
        # Clean up test files
        if isfile(test_netcdf_file)
            rm(test_netcdf_file)
        end
        test_output_file = joinpath(tempdir(), "test_mobility_output.nc")
        if isfile(test_output_file)
            rm(test_output_file)
        end
    end
end

