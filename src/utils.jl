## ----------------------------------------------------------------------------
## EPIDEMIC PARAMS RELATED FUNCTIONS
## ----------------------------------------------------------------------------

"""
Epidemic_Params
  Struct that contains the parameters related with the epidemic parameters and compartmental evolution.

# Fields
  All the parameters contained in this structure are probabilities ranged between 0 and 1.

# Epidemic parameters
  - `βᴵ::Array{Float64, 1}`: Infectivity of infected.
  - `βᴬ::Array{Float64, 1}`: Infectivity of asymptomatic.
  - `ηᵍ::Array{Float64, 1}`: Exposed rate for each strata.
  - `αᵍ::Array{Float64, 1}`: Asymptomatic infectious rate for each strata.
  - `μᵍ::Array{Float64, 1}`: Infectious rate for each strata.
  - `θᵍ::Array{Float64, 2}`: Direct death probability for each strata.
  - `γᵍ::Array{Float64, 2}`: ICU probability for each strata.
  - `ζᵍ::Array{Float64, 1}`: Pre-deceased rate for each strata.
  - `λᵍ::Array{Float64, 1}`: Pre-hospitalized in ICU rate for each strata.
  - `ωᵍ::Array{Float64, 2}`: Fatality probability in ICU for each strata.
  - `ψᵍ::Array{Float64, 1}`: Death rate in iCU for each strata.
  - `χᵍ::Array{Float64, 1}`: ICU discharge rate for each strata.
  - `Λ::Float64`: Probability of losing the vaccine-acquired immunity
  - `Γ::Float64`: Probability of losing the disease-acquired immunity
  - `T::Int64`: Number of epidemic timesteps.
  - `V::Int64`: Number of stages in the vaccination process
  - `NumComps::Int64`: Number of epidemiological compartments
  - `rᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing infections
  - `kᵥ::Array{Float64, 1}`: Vaccine efficacy in preventing tranmission

# Compartmental evolution
  - `ρˢᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of suceptible individuals for each
    strata, patch and vaccination status.
  - `ρᴱᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of exposed individuals for each
    strata, patch and vaccination status.
  - `ρᴬᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of asymptomatic individuals for
    each strata, patch and vaccination status.
  - `ρᴵᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of symptomatic individuals for each
    strata, patch and vaccination status.
  - `ρᴾᴴᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of pre-hospitalized to ICU
    individuals for each strata, patch and vaccination status.
  - `ρᴾᴰᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of pre-deceased individuals for each strata, patch 
    and vaccination status.
  - `ρᴴᴿᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of hospitalized in ICU patients who
    will recover for each strata, patch and vaccination status.
  - `ρᴴᴰᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of hospitalized in ICU patients who
    will not recover for each strata, patch and vaccination status.
  - `ρᴰᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of deceased individuals for each
    strata, patch and vaccination status.
  - `ρᴿᵍᵥ::Array{Float64, 4}`: Matrix of size ``G \\times M \\times T \\times V`` containing
    infomation about the evolution of fraction of recovered individuals for each
    strata, patch and vaccination status.

# Auxiliary
  - `CHᵢᵍᵥ::Array{Float64, 3}`: Fraction of securely confined individuals for each
    strata and patch.
  - `Qᵢᵍ::Array{Float64, 3}`: Suceptible contacts available for each strata on a
    given patch.

"""
struct Epidemic_Params
    #Epidemic parameters
    βᴵ::Array{Float64, 1}
    βᴬ::Array{Float64, 1}
    ηᵍ::Array{Float64, 1}
    αᵍ::Array{Float64, 1}
    μᵍ::Array{Float64, 1}
    θᵍ::Array{Float64, 2}
    γᵍ::Array{Float64, 2}
    ζᵍ::Array{Float64, 1}
    λᵍ::Array{Float64, 1}
    ωᵍ::Array{Float64, 2}
    ψᵍ::Array{Float64, 1}
    χᵍ::Array{Float64, 1}
    Λ::Float64
    Γ::Float64
    
    # Epidemic parameter for vaccinations
    rᵥ::Array{Float64, 1}
    kᵥ::Array{Float64, 1}

    # The total number of agents, patches, time-steps, 
    # vaccination-states, and compartments
    G::Int64
    M::Int64
    T::Int64
    V::Int64
    NumComps::Int64

    CompLabels::Array{String, 1}
    VaccLabels::Array{String, 1}

    # Compartments evolution
    ρˢᵍᵥ::Array{Float64, 4}
    ρᴱᵍᵥ::Array{Float64, 4}
    ρᴬᵍᵥ::Array{Float64, 4}
    ρᴵᵍᵥ::Array{Float64, 4}
    ρᴾᴴᵍᵥ::Array{Float64, 4}
    ρᴾᴰᵍᵥ::Array{Float64, 4}
    ρᴴᴿᵍᵥ::Array{Float64, 4}
    ρᴴᴰᵍᵥ::Array{Float64, 4}
    ρᴰᵍᵥ::Array{Float64, 4}
    ρᴿᵍᵥ::Array{Float64, 4}
    
    # Fraction of securely confined individuals for each strata and patch.
    CHᵢᵍᵥ::Array{Float64, 4}
    
    # R_t related arrays
    Qᵢᵍ::Array{Float64, 3}
end


"""
    Epidemic_Params(βᴵ::Float64,
                    βᴬ::Float64,
                    ηᵍ::Array{Float64, 1},
                    αᵍ::Array{Float64, 1},
                    μᵍ::Array{Float64, 1},
                    θᵍ::Array{Float64, 2},
                    γᵍ::Array{Float64, 2},
                    ζᵍ::Array{Float64, 1},
                    λᵍ::Array{Float64, 1},
                    ωᵍ::Array{Float64, 2},
                    ψᵍ::Array{Float64, 1},
                    χᵍ::Array{Float64, 1},
                    Λ::Float64,
                    Γ::Float64,                         
                    rᵥ::Array{Float64, 1},
                    kᵥ::Array{Float64, 1},
                    G::Int64,
                    M::Int64,
                    T::Int64,
                    V::Int64)

# Arguments
  - `βᴵ::Float64`: Infectivity of infected for each vaccination status.
  - `βᴬ::Float64`: Infectivity of asymptomatic for each vaccination status.
  - `ηᵍ::Array{Float64, 1}`: Vector of size ``G`` with exposed rates for each
    strata.
  - `αᵍ::Array{Float64, 1}`: Vector of size ``G`` with asymptomatic infectious
    rates for each strata.
  - `μᵍ::Array{Float64, 1}`: Vector of size ``G`` with infectious rates for each
    strata.
  - `θᵍ::Array{Float64, 2}`: Matrix of size ``GxV`` with direct death probabilities
    for each strata.
  - `γᵍ::Array{Float64, 2}`: Matrix of size ``GxV`` with ICU probabilities for each
    strata.
  - `ζᵍ::Array{Float64, 1}`: Vector of size ``G`` with pre-deceased rates for
    each strata.
  - `λᵍ::Array{Float64, 1}`: Vector of size ``G`` with pre-hospitalized in ICU
    rates for each strata.
  - `ωᵍ::Array{Float64, 2}`: Matrix of size ``GxV`` with fatality probabilities in
    ICU for each strata.
  - `ψᵍ::Array{Float64, 1}`: Vector of size ``G`` with death rates for each
    strata.
  - `χᵍ::Array{Float64, 1}`: Vector of size ``G`` with ICU discharge rates for
    each strata.
  - `Λ::Float64`: Probability of losing the vaccine-acquired immunity
  - `Γ::Float64`: Probability of losing the disease-acquired immunity
  - `rᵥ::Array{Float64, 1}`: Vector of size ``V``. Vaccine efficacy in preventing infections 
    for each vaccination status.
  - `kᵥ::Array{Float64, 1}`: Vector of size ``V``. Vaccine efficacy in preventing tranmission 
    for each vaccination status.
  - `G::Int64`: Number of strata.
  - `M::Int64`: Number of patches.
  - `T::Int64`: Number of epidemic timesteps.
  - `V::Int64`: Number of stages in the vaccination process.

# Return

Struct that contains the parameters related with the epidemic parameters and
compartmental evolution.
"""
function Epidemic_Params(βᴵ::Float64,
                         βᴬ::Float64,
                         ηᵍ::Array{Float64, 1},
                         αᵍ::Array{Float64, 1},
                         μᵍ::Array{Float64, 1},
                         θᵍ::Array{Float64, 2},
                         γᵍ::Array{Float64, 2},
                         ζᵍ::Array{Float64, 1},
                         λᵍ::Array{Float64, 1},
                         ωᵍ::Array{Float64, 2},
                         ψᵍ::Array{Float64, 1},
                         χᵍ::Array{Float64, 1},
                         Λ::Float64,
                         Γ::Float64,                         
                         rᵥ::Array{Float64, 1},
                         kᵥ::Array{Float64, 1},
                         G::Int64,
                         M::Int64,
                         T::Int64)

    NumComps = 10
    V = 3
    CompLabels = ["S", "E", "A", "I", "PH", "PD", "HR", "HD", "R", "D"]
    VaccLabels = ["NV", "V", "PV"]
    # Allocate memory for simulations
    ρˢᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴱᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴬᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴵᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴾᴴᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴾᴰᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴴᴿᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴴᴰᵍᵥ = zeros(Float64, G, M, T, V)
    ρᴰᵍᵥ  = zeros(Float64, G, M, T, V)
    ρᴿᵍᵥ  = zeros(Float64, G, M, T, V)
    CHᵢᵍᵥ = zeros(Float64, G, M, T, V)
    Qᵢᵍ   = zeros(Float64, G, M, T)
    
    return Epidemic_Params([βᴵ], [βᴬ], copy(ηᵍ), copy(αᵍ), copy(μᵍ),
                           copy(θᵍ), copy(γᵍ), copy(ζᵍ), copy(λᵍ), copy(ωᵍ),
                           copy(ψᵍ), copy(χᵍ), copy(Λ), copy(Γ), copy(rᵥ), copy(kᵥ), 
                           G, M, T, V, NumComps, CompLabels, VaccLabels, 
                           ρˢᵍᵥ, ρᴱᵍᵥ, ρᴬᵍᵥ, ρᴵᵍᵥ, ρᴾᴴᵍᵥ, ρᴾᴰᵍᵥ, ρᴴᴿᵍᵥ, ρᴴᴰᵍᵥ, ρᴰᵍᵥ, 
                           ρᴿᵍᵥ, CHᵢᵍᵥ, Qᵢᵍ)
end


function update_epidemic_params!(epi_params::Epidemic_Params, data_dict::Dict{String, Any})
  # has_key(key) = haskey(data_dict, key) && !(typeof(getfield(person, Symbol(key))) <: Union{Missing, Nothing})
  epi_params.βᴵ = has_key("βᴵ") ? data_dict["βᴵ"] : epi_params.βᴵ
  epi_params.βᴬ = has_key("βᴬ") ? data_dict["βᴬ"] : epi_params.βᴬ
  epi_params.ηᵍ = has_key("ηᵍ") ? data_dict["ηᵍ"] : epi_params.ηᵍ
  epi_params.αᵍ = has_key("αᵍ") ? data_dict["αᵍ"] : epi_params.αᵍ
  epi_params.μᵍ = has_key("μᵍ") ? data_dict["μᵍ"] : epi_params.μᵍ
  epi_params.θᵍ = has_key("θᵍ") ? data_dict["θᵍ"] : epi_params.θᵍ
  epi_params.γᵍ = has_key("γᵍ") ? data_dict["γᵍ"] : epi_params.γᵍ
  epi_params.ζᵍ = has_key("ζᵍ") ? data_dict["ζᵍ"] : epi_params.ζᵍ
  epi_params.λᵍ = has_key("λᵍ") ? data_dict["λᵍ"] : epi_params.λᵍ
  epi_params.ωᵍ = has_key("ωᵍ") ? data_dict["ωᵍ"] : epi_params.ωᵍ
  epi_params.ψᵍ = has_key("ψᵍ") ? data_dict["ψᵍ"] : epi_params.ψᵍ
  epi_params.χᵍ = has_key("χᵍ") ? data_dict["χᵍ"] : epi_params.χᵍ
  epi_params.Λ  = has_key("Λ")  ? data_dict["Λ"] : epi_params.Λ
  epi_params.Γ  = has_key("Γ")  ? data_dict["Γ"] : epi_params.Γ
  epi_params.rᵥ = has_key("rᵥ") ? data_dict["rᵥ"] : epi_params.rᵥ
  epi_params.kᵥ = has_key("kᵥ") ? data_dict["kᵥ"] : epi_params.kᵥ
end


"""
    reset_epidemic_compartments!(epi_params::Epidemic_Params)
Reset the ρ's to reuse the structure and avoid additional allocations.

# Arguments
- `epi_params::Epidemic_Params`: structure to reset.
"""
function reset_epidemic_compartments!(epi_params::Epidemic_Params)
    epi_params.ρˢᵍᵥ .= 0.
    epi_params.ρˢᵍᵥ[:, :, :, 1] .= 1.
    epi_params.ρᴱᵍᵥ  .= 0.
    epi_params.ρᴬᵍᵥ  .= 0.
    epi_params.ρᴵᵍᵥ  .= 0.
    epi_params.ρᴾᴴᵍᵥ .= 0.
    epi_params.ρᴾᴰᵍᵥ .= 0.
    epi_params.ρᴴᴿᵍᵥ .= 0.
    epi_params.ρᴴᴰᵍᵥ .= 0.
    epi_params.ρᴰᵍᵥ  .= 0.
    epi_params.ρᴿᵍᵥ  .= 0.
    epi_params.CHᵢᵍᵥ .= 0.

    # Rt structures
    epi_params.Qᵢᵍ .= 0.

    nothing
end

### ----------------------------------------------------------------------------
### PATCH AND POPULATION RELATED FUNCTIONS
### ----------------------------------------------------------------------------

"""
Population_Params

Struct that contains the parameters related with geographical, population and
mobility data.

# Fields
  - `G::Int64`: Number of population strata.
  - `M::Int64`: Number of patches.
  - `nᵢ:Array{Float64, 1}`: Population of each patch.
  - `nᵢᵍ:Array{Float64, 2}`: Population of each strata on each patch.
  - `nᵢ_eff:Array{Float64, 1}`: Effective population on each patch after taking
    into account mobility.
  - `nᵢᵍ_eff:Array{Float64, 1}`: Effective population of each strata on each patch
    after taking into account mobility.
  - `N::Int64`: Total population.
  - `Nᵍ::Array{Int64, 1}`: Total population of each strata.
  - `pᵍ::Array{Float64, 1}`: Vector with the degree of mobility of each strata.
  - `pᵍ_eff::Array{Float64, 1}`: Vector with the current degree of mobility of
    each strata.
  - `edgelist::Array{Int64, 2}`: Matrix with the directed edgelist between
    patches, where ``L`` is the number of edges.
  - `Rᵢⱼ::Array{Float64, 1}`: Vector with the transition probabilities for each
    edge in the edgelist.
  - `mobilityᵍ::Array{Float64, 1}`: Effective mobility of each strata between
    patches (used to optimize performance).
  - `kᵍ::Array{Float64, 1}`: Average number of contacts of each strata.
  - `kᵍ_h::Array{Float64, 1}`: Average number of contacts at home of each strata.
  - `kᵍ_w::Array{Float64, 1}`: Average number of contacts at work of each strata.
  - `kᵍ_eff::Array{Float64, 1}`: Current average number of contacts of each
    strata.
  - `C::Array{Float64, 2}`: Matrix with the probability of contact between
    different stratifications.
  - `zᵍ::Array{Float64, 1}`: Nomalization factor each strata.
  - `normᵍ::Array{Float64, 2}`: Normalization of each strata (used to optimize
    performance).
  - `sᵢ::Array{Float64, 1}`: Surface of each patch.
  - `ξ::Float64`: Densisty factor.
  - `σ::Float64`: Average household size.

# Constructor

Use outer constructor [`Population_Params`](#MMCAcovid19.Population_Params-Tuple{Int64,Int64,Array{Float64,2},Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,2},Array{Float64,1},Array{Int64,2},Array{Float64,1},Array{Float64,1},Float64,Float64})
for a proper initialization of this struct.
"""
struct Population_Params
    G::Int64
    M::Int64
    nᵢ::Array{Float64, 1}
    nᵢᵍ::Array{Float64, 2}
    nᵢ_eff::Array{Float64, 1}
    nᵢᵍ_eff::Array{Float64, 2}
    N::Int64
    Nᵍ::Array{Int64, 1}
    pᵍ::Array{Float64, 1}
    pᵍ_eff::Array{Float64, 1}
    edgelist::Array{Int64, 2}
    Rᵢⱼ::Array{Float64, 1}
    mobilityᵍ::Array{Float64, 2}
    kᵍ::Array{Float64, 1}
    kᵍ_h::Array{Float64, 1}
    kᵍ_w::Array{Float64, 1}
    kᵍ_eff::Array{Float64, 1}
    C::Array{Float64, 2}
    zᵍ::Array{Float64, 1}
    normᵍ::Array{Float64, 2}
    sᵢ::Array{Float64, 1}
    ξ::Float64
    σ::Float64
end


"""
    Population_Params(G::Int64,
                      M::Int64,
                      nᵢᵍ::Array{Float64, 2},
                      kᵍ::Array{Float64, 1},
                      kᵍ_h::Array{Float64, 1},
                      kᵍ_w::Array{Float64, 1},
                      C::Array{Float64, 2},
                      pᵍ::Array{Float64, 1},
                      edgelist::Array{Int64, 2},
                      Rᵢⱼ::Array{Float64, 1},
                      sᵢ::Array{Float64, 1},
                      ξ::Float64,
                      σ::Float64)

Constructor of the struct [`Population_Params`](@ref).

# Arguments
  - `G::Int64`: Number of population strata.
  - `M::Int64`: Number of patches.
  - `nᵢᵍ::Array{Float64, 2}`: Matrix of size ``G \\times M`` with the population
    at each strata and patch.
  - `kᵍ::Array{Float64, 1}`: Vector of size ``G`` with the average number of
    contacts of each strata.
  - `kᵍ_h::Array{Float64, 1}`: Vector of size ``G`` with the average number of
    contacts at home of each strata.
  - `kᵍ_w::Array{Float64, 1}`: Vector of size ``G`` with the average number of
    contacts at work of each strata.
  - `C::Array{Float64, 2}`: Matrix of size ``G \\times G`` with the probability of
    contact between different stratifications.
  - `pᵍ::Array{Float64, 1}`: Vector of size ``G`` with the degree of mobility of
    each strata ranged between 0 and 1.
  - `edgelist::Array{Int64, 2}`: Matrix of size ``L \\times 2`` containing the
    directed edgelist between patches, where ``L`` is the number of edges. The IDs
    of the patches have to go from ``1`` to ``M``.
  - `Rᵢⱼ::Array{Float64, 1}`: Vector of size ``L`` containing the transition
    probabilities for each edge in the edgelist.
  - `sᵢ::Array{Float64, 1}`: Vector of size `M` with the surface of each patch.
  - `ξ::Float64`: Density factor.
  - `σ::Float64`: Average household size.

# Return
  - Struct that contains the parameters related with geographical, population and mobility data.
"""
function Population_Params(G::Int64,
                           M::Int64,
                           nᵢᵍ::Array{Float64, 2},
                           kᵍ::Array{Float64, 1},
                           kᵍ_h::Array{Float64, 1},
                           kᵍ_w::Array{Float64, 1},
                           C::Array{Float64, 2},
                           pᵍ::Array{Float64, 1},
                           edgelist::Array{Int64, 2},
                           Rᵢⱼ::Array{Float64, 1},
                           sᵢ::Array{Float64, 1},
                           ξ::Float64,
                           σ::Float64)

    edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)

    # Aggregate population count
    Nᵍ = round.(Int, sum(nᵢᵍ, dims = 2)[:, 1])
    N = sum(Nᵍ)
    nᵢ = sum(nᵢᵍ, dims = 1)[1, :]
    mobilityᵍ = zeros(Float64, G, length(Rᵢⱼ))

    # Init. effective population
    nᵢ_eff = zeros(Float64, M)
    nᵢᵍ_eff = (1 .- pᵍ) .* nᵢᵍ

    # Init. normalization vector
    zᵍ = zeros(Float64, G)
    normᵍ = zeros(Float64, G, M)

    # Compute effective population
    compute_effective_population!(nᵢᵍ_eff, nᵢ_eff, nᵢᵍ, Nᵍ, mobilityᵍ, kᵍ, zᵍ,
                                  normᵍ, ξ, pᵍ, sᵢ, edgelist, Rᵢⱼ, M, G)

    return Population_Params(G, M, nᵢ, copy(nᵢᵍ), nᵢ_eff, nᵢᵍ_eff, N, Nᵍ,
                             copy(pᵍ), copy(pᵍ), edgelist, Rᵢⱼ, mobilityᵍ,
                             copy(kᵍ), copy(kᵍ_h), copy(kᵍ_w), copy(kᵍ), copy(C),
                             zᵍ, normᵍ, copy(sᵢ), ξ, σ)
end



"""
    reset_params!(epi_params::Epidemic_Params,
                  population::Population_Params)

Reset the epidemic and population parameters to reuse the structure and avoid
additional data allocations.

# Arguments
  - `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
    and the epidemic spreading information.
  - `population::Population_Params`: Structure that contains all the parameters
  related with the population.
"""
function reset_params!(epi_params::Epidemic_Params,
                       population::Population_Params)

    reset_epidemic_params!(epi_params)
    population.pᵍ_eff .= population.pᵍ
    population.kᵍ_eff .= population.kᵍ
    update_population_params!(population)

end

### ---------------------------------------------------------------
### CONTFINAMENT AND VACCINATION PARAMETERS
### ---------------------------------------------------------------


"""
Struct that contains the parameters related with containement
measures

#Parameters
  - κ₀s:Array{Float64, 1}: Array of level of confinement
  - ϕs:Array{Float64, 1}: Array of permeabilities of confined households
  - δs:Array{Float64, 1}: Array of social distancing measures
  - tᶜs:Array{Int64, 1}: Timesteps when the containment measures will be applied
"""

struct NPI_Params
    κ₀s::Array{Float64, 1}
    ϕs::Array{Float64, 1}
    δs::Array{Float64, 1}
    tᶜs::Array{Int64, 1}
end



"""
VACCINATION PARAMETERS

Struct that contains the parameters related with vaccination 
strageties

#Parameters
  - ϵᵍ: total vaccinations per age strata
  - tᵛs: 
  - ϵᵍs:

"""




### ---------------------------------------------------------------
### FUNCTION TO INITIALIZE PARAMETER STRUCTS FROM DICTS
### ---------------------------------------------------------------

function init_pop_param_struct(G::Int64, M::Int64,
                               G_coords::Array{String, 1},
                               pop_params_dict::Dict, 
                               metapop_df::DataFrame,
                               network_df::DataFrame)

    # Subpopulations' patch surface
    sᵢ = metapop_df[:, "area"]
    # Subpopulation by age strata
    nᵢᵍ = copy(transpose(Array{Float64,2}(metapop_df[:, G_coords])))
    # Age Contact Matrix
    C = Float64.(mapreduce(permutedims, vcat, pop_params_dict["C"]))
    # Average number of contacts per strata
    kᵍ = Float64.(pop_params_dict["kᵍ"])
    # Average number of contacts at home per strata
    kᵍ_h = Float64.(pop_params_dict["kᵍ_h"])
    # Average number of contacts at work per strata
    kᵍ_w = Float64.(pop_params_dict["kᵍ_w"])
    # Degree of mobility per strata
    pᵍ = Float64.(pop_params_dict["pᵍ"])
    # Density factor
    ξ = pop_params_dict["σ"]
    # Average household size
    σ = pop_params_dict["σ"]

    edgelist = Array{Int64, 2}(network_df[:, 1:2])
    Rᵢⱼ      = copy(network_df[:, 3])
    edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)
    pop_params    = Population_Params(G, M, nᵢᵍ, kᵍ, kᵍ_h, kᵍ_w, C, pᵍ, edgelist, Rᵢⱼ, sᵢ, ξ, σ)

    return pop_params
end

function init_epi_parameters_struct(G::Int64, M::Int64, T::Int64,
                                    G_coords::Array{String, 1}, 
                                    epi_params_dict::Dict)

    # Scaling of the asymptomatic infectivity
    scale_β = epi_params_dict["scale_β"]
    # Infectivity of Symptomatic
    βᴵ = epi_params_dict["βᴵ"]
    # Infectivity of Asymptomatic
    if haskey(epi_params_dict, "βᴬ")
        βᴬ = epi_params_dict["βᴬ"]
    elseif haskey(epi_params_dict, "scale_β")
        βᴬ = scale_β * βᴵ
    else
        @error "Either βᴬ or scale_β should be provided"
    end
    # Exposed rate
    ηᵍ = Float64.(epi_params_dict["ηᵍ"])
    # Asymptomatic rate
    αᵍ = Float64.(epi_params_dict["αᵍ"])
    # Infectious rate
    μᵍ = Float64.(epi_params_dict["μᵍ"])

    # Waning immunity rate 
    Λ = epi_params_dict["Λ"] 
    # Reinfection rate
    Γ = epi_params_dict["Γ"] 


    ## EPIDEMIC PARAMETERS TRANSITION RATES VACCINATION

    # Direct death probability
    θᵍ = Float64.(reduce(hcat, [epi_params_dict["θᵍ"], epi_params_dict["θᵍ"] * epi_params_dict["risk_reduction_dd"], 
                                epi_params_dict["θᵍ"] * epi_params_dict["risk_reduction_dd"]]) )
    # Hospitalization probability
    γᵍ = Float64.(reduce(hcat, [epi_params_dict["γᵍ"], epi_params_dict["γᵍ"] * epi_params_dict["risk_reduction_h"],
                                epi_params_dict["γᵍ"] * epi_params_dict["risk_reduction_h"]]) )
    # Fatality probability in ICU
    ωᵍ = Float64.(reduce(hcat, [epi_params_dict["ωᵍ"], epi_params_dict["ωᵍ"] * epi_params_dict["risk_reduction_d"], 
                                epi_params_dict["ωᵍ"] * epi_params_dict["risk_reduction_d"]]) )
    # Pre-deceased rate
    ζᵍ = Float64.(epi_params_dict["ζᵍ"])
    # Pre-hospitalized in ICU rate
    λᵍ = Float64.(epi_params_dict["λᵍ"])
    # Death rate in ICU
    ψᵍ = Float64.(epi_params_dict["ψᵍ"])
    # ICU discharge rate
    χᵍ = Float64.(epi_params_dict["χᵍ"])
    # Relative risk reduction of the probability of infection
    rᵥ = Float64.(epi_params_dict["rᵥ"])
    # Relative risk reduction of the probability of transmission
    kᵥ = Float64.(epi_params_dict["kᵥ"])

    return Epidemic_Params(βᴵ, βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ, ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ, Λ, Γ, rᵥ, kᵥ, G, M, T)
end

function init_NPI_parameters_struct(data_path::String, npi_params_dict::Dict, kappa0_filename::String, first_day::Date)
    κ₀_df = nothing
    if !isnothing(kappa0_filename)
        kappa0_filename = joinpath(data_path, kappa0_filename)
        @info "- Loading κ₀ time series from $(kappa0_filename)"
        κ₀_df = CSV.read(kappa0_filename, DataFrame);
    end

    return init_NPI_parameters_struct(κ₀_df, npi_params_dict, first_day)
end

function init_NPI_parameters_struct(κ₀_df::Union{DataFrame, Nothing}, npi_params_dict::Dict, first_day::Date)
    if isnothing(κ₀_df)
        # Handle the case when κ₀_df is nothing
        # For example, you could use default values or read from npi_params_dict
        tᶜs = npi_params_dict["tᶜs"]
        κ₀s = npi_params_dict["κ₀s"]
        ϕs = npi_params_dict["ϕs"]
        δs = npi_params_dict["δs"]
    else
        # Existing logic for when κ₀_df is a DataFrame
        @info "- Synchronizing to dates"
        κ₀_df.time = map(x -> (x .- first_day).value + 1, κ₀_df.date)
        # Timesteps when the containment measures will be applied
        tᶜs = κ₀_df.time[:]
        # Array of level of confinement
        κ₀s = κ₀_df.reduction[:]
        # Array of premeabilities of confined households
        
        ϕs_aux = Float64.(npi_params_dict["ϕs"])
        δs_aux = Float64.(npi_params_dict["δs"])

        # Supposing ϕs and δs are constant, while the confinement measures are applied
        ϕs = fill(ϕs_aux[1], length(tᶜs))
        # Array of social distancing measures
        δs = fill(δs_aux[1], length(tᶜs))
    end

    return NPI_Params(κ₀s, ϕs, δs, tᶜs)
end


### ---------------------------------------------------------------
### FUNCTION TO SET INITIAL CONDITIONS
### ---------------------------------------------------------------

function set_compartments!(epi_params, population, 
                            initial_compartments; normalize=true)
    V = epi_params.V
    # Index of the initial condition
    t₀ = 1
    if normalize
        for i in 1:V
            epi_params.ρˢᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 1] ./ population.nᵢᵍ
            epi_params.ρᴱᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 2] ./ population.nᵢᵍ
            epi_params.ρᴬᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 3] ./ population.nᵢᵍ
            epi_params.ρᴵᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 4] ./ population.nᵢᵍ
            epi_params.ρᴾᴴᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 5] ./ population.nᵢᵍ
            epi_params.ρᴾᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 6] ./ population.nᵢᵍ
            epi_params.ρᴴᴿᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 7] ./ population.nᵢᵍ
            epi_params.ρᴴᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 8] ./ population.nᵢᵍ
            epi_params.ρᴿᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 9] ./ population.nᵢᵍ
            epi_params.ρᴰᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 10] ./ population.nᵢᵍ
        end
    else
        for i in 1:V
            epi_params.ρˢᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 1] 
            epi_params.ρᴱᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 2] 
            epi_params.ρᴬᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 3] 
            epi_params.ρᴵᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 4] 
            epi_params.ρᴾᴴᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 5] 
            epi_params.ρᴾᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 6] 
            epi_params.ρᴴᴿᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 7] 
            epi_params.ρᴴᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 8] 
            epi_params.ρᴿᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 9] 
            epi_params.ρᴰᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 10]
        end
    end

    epi_params.ρˢᵍᵥ[isnan.(epi_params.ρˢᵍᵥ)]   .= 0
    epi_params.ρᴱᵍᵥ[isnan.(epi_params.ρᴱᵍᵥ)]   .= 0
    epi_params.ρᴬᵍᵥ[isnan.(epi_params.ρᴬᵍᵥ)]   .= 0
    epi_params.ρᴵᵍᵥ[isnan.(epi_params.ρᴵᵍᵥ)]   .= 0
    epi_params.ρᴾᴴᵍᵥ[isnan.(epi_params.ρᴾᴴᵍᵥ)] .= 0
    epi_params.ρᴾᴰᵍᵥ[isnan.(epi_params.ρᴾᴰᵍᵥ)] .= 0
    epi_params.ρᴴᴿᵍᵥ[isnan.(epi_params.ρᴴᴿᵍᵥ)] .= 0
    epi_params.ρᴴᴰᵍᵥ[isnan.(epi_params.ρᴴᴰᵍᵥ)] .= 0
    epi_params.ρᴿᵍᵥ[isnan.(epi_params.ρᴿᵍᵥ)]   .= 0
    epi_params.ρᴰᵍᵥ[isnan.(epi_params.ρᴰᵍᵥ)]   .= 0
end


"""
set_initial_conditions!(epi_params::Epidemic_Params,
     population::Population_Params,
     Sᵛ₀::Array{Float64, 2},
     E₀::Array{Float64, 2},
     A₀::Array{Float64, 2},
     I₀::Array{Float64, 2},
     H₀::Array{Float64, 2},
     R₀::Array{Float64, 2})

Set the initial number of individuals on a population. They can be
introduced as susceptible vaccinated (Sᵛ₀), exposed unvaccinated (E₀), asymptomatic 
unvaccinated (A₀), symptomatic unvaccinated (I₀), hospitalized unvaccinated (H₀) or recovered 
unvaccinated (R₀) individuals.

# Arguments
- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
related with the population.
- `Sᵛ₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of vaccinated susceptible individuals of each strata on each patch.
- `E₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of unvaccinated exposed individuals of each strata on each patch.
- `A₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of unvaccinated asymptomatic infected individuals of each strata on each patch.
- `I₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of unvaccinated symptomatic infected individuals of each strata on each patch.
- `H₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of unvaccinated pre-hospitalized individuals of each strata on each patch.
- `R₀::Array{Float64, 2}`: Matrix of size ``G \\times M`` containing the number
of unvaccinated recovered individuals of each strata on each patch.
"""
function set_initial_conditions!(epi_params::Epidemic_Params,
                               population::Population_Params,
                               Sᵛ₀::Array{Float64, 2},
                               E₀::Array{Float64, 2},
                               A₀::Array{Float64, 2},
                               I₀::Array{Float64, 2},
                               H₀::Array{Float64, 2},
                               R₀::Array{Float64, 2})

    t₀ = 1
    
    # Initial susceptible population
    @. epi_params.ρˢᵍᵥ[:, :, t₀, 2] = Sᵛ₀ / population.nᵢᵍ

    # Initial exposed population
    @. epi_params.ρᴱᵍᵥ[:, :, t₀, 1] = E₀ / population.nᵢᵍ 

    # Initial asymptomatic population
    @. epi_params.ρᴬᵍᵥ[:, :, t₀, 1] = A₀ / population.nᵢᵍ 

    # Initial sypmtomatic population
    @. epi_params.ρᴵᵍᵥ[:, :, t₀, 1] = I₀ / population.nᵢᵍ 
    
    # Initial hospitalized population
    @. epi_params.ρᴾᴴᵍᵥ[:, :, t₀, 1] = H₀ / population.nᵢᵍ 
    
    # Initial recovered population
    @. epi_params.ρᴿᵍᵥ[:, :, t₀, 1] = R₀ / population.nᵢᵍ

    # Control over division by zero
    epi_params.ρˢᵍᵥ[isnan.(epi_params.ρˢᵍᵥ)] .= 0
    epi_params.ρᴱᵍᵥ[isnan.(epi_params.ρᴱᵍᵥ)] .= 0
    epi_params.ρᴬᵍᵥ[isnan.(epi_params.ρᴬᵍᵥ)] .= 0
    epi_params.ρᴵᵍᵥ[isnan.(epi_params.ρᴵᵍᵥ)] .= 0
    epi_params.ρᴾᴴᵍᵥ[isnan.(epi_params.ρᴾᴴᵍᵥ)] .= 0
    epi_params.ρᴿᵍᵥ[isnan.(epi_params.ρᴿᵍᵥ)] .= 0

    # Update the fraction of suceptible individual
    @. epi_params.ρˢᵍᵥ[:, :, t₀, 1] = 1 - (epi_params.ρˢᵍᵥ[:, :, t₀, 2] +
                                       epi_params.ρᴱᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴬᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴵᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴾᴴᵍᵥ[:, :, t₀, 1] +
                                       epi_params.ρᴿᵍᵥ[:, :, t₀, 1])
    
    @. epi_params.ρˢᵍᵥ[abs(epi_params.ρˢᵍᵥ) < 1e-12]  = 0
    
    nothing
end




### ----------------------------------------------------------------------------
### EXTRA FUNCTIONS
### ----------------------------------------------------------------------------

"""
    correct_self_loops(edgelist::Array{Int64, 2},
                       Rᵢⱼ::Array{Float64, 1},
                       M::Int64)

Repair weighted matrices to ensure that all nodes are represented adding an
additional self-loop.

# Arguments
- `edgelist::Array{Int64, 2}`: Matrix containing the directed edgelist between
  patches. The IDs of the patches have to go from 1 to M.
- `Rᵢⱼ::Array{Float64, 1}`: Vector containing the transition probabilities for
  each edge in the edgelist.
- `M::Int64`: Number of patches.

# Return

Returns the corrected list of edges and transition probabilities.
"""
function correct_self_loops(edgelist::Array{Int64, 2},
                            Rᵢⱼ::Array{Float64, 1},
                            M::Int64)

    self_loops = falses(M)
    k_in = zeros(Int64, M)
    k_out = zeros(Int64, M)

    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        k_in[j] += 1
        k_out[i] += 1

        if i == j
            self_loops[i] = true
        end
    end

    to_add = collect(1:M)[(k_out .== 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, ones(Float64, length(to_add)))

    to_add = collect(1:M)[(k_out .!= 0) .& .!(self_loops)]

    new_edges = zeros(Int64, length(to_add), 2)
    new_edges[:,1] = to_add
    new_edges[:,2] = to_add

    edgelist = vcat(edgelist, new_edges)
    Rᵢⱼ = vcat(Rᵢⱼ, zeros(Float64, length(to_add)))

    return (edgelist, Rᵢⱼ)
end

"""
    make_edls(Rᵢⱼ) 

Turns a origin destination matrix into an edgelist

# Arguments
- `Rᵢⱼ::Array{Float64, 2}`: Matrix containing the flows of people from and to each node
  of the networks

# Return

Returns an edgelist and a list of flows in those edges.
"""

function make_edls(Rᵢⱼ) 
    Rᵢⱼ = replace(Rᵢⱼ, NaN => 0)
    M = size(Rᵢⱼ)[1]
    edgl = ones( sum( Rᵢⱼ .!= 0 ), 2 )
    count = 1
    a = ones( sum(Rᵢⱼ .!= 0) )
    
    for i in 1:M
        @simd for j in 1:M
            if Rᵢⱼ[i, j] != 0
                edgl[count, 1] = i
                edgl[count, 2] = j
                a[count] = Rᵢⱼ[i, j]
                count += 1
            end
        end
    end  
    
    edgl = convert(Matrix{Int64}, edgl)
    
    return edgl, a
end
