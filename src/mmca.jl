"""
    update_prob!(Pᵢᵍᵥ::Array{Float64, 3},
                 Sᵢᵍᵥ::Array{Float64, 3},
                 τᵢᵍᵥ::Array{Float64, 3},
                 epi_params::Epidemic_Params,
                 population::Population_Params,
                 κ₀::Float64,
                 ϕ::Float64,
                 δ::Float64,
                 ϵᵍ::Array{Float64, 1},
                 t::Int64,
                 tᶜ::Int64,
                 tᵛ::Int64)

Updates the probabilities of the model using the equations described in the
paper.
"""
function update_prob!(Pᵢᵍᵥ::Array{Float64, 3},
                      Sᵢᵍᵥ::Array{Float64, 3},
                      τᵢᵍᵥ::Array{Float64, 3},
                      epi_params::Epidemic_Params,
                      population::Population_Params,
                      κ₀::Float64,
                      ϕ::Float64,
                      δ::Float64,
                      ϵᵍ::Array{Float64, 1},
                      t::Int64,
                      tᶜ::Int64, 
                      tᵛ::Int64)

    # Shortcuts to parameters
    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    θᵍᵥ = epi_params.θᵍᵥ
    γᵍᵥ = epi_params.γᵍᵥ
    ζᵍ = epi_params.ζᵍ
    λᵍ = epi_params.λᵍ
    ωᵍᵥ = epi_params.ωᵍᵥ
    ψᵍ = epi_params.ψᵍ
    χᵍ = epi_params.χᵍ
    rᵥ = epi_params.rᵥ
    kᵥ = epi_params.kᵥ
    Γ = epi_params.Γ
    Λ = epi_params.Λ
    
    # Shortcut to variables
    ρˢᵍᵥ = epi_params.ρˢᵍᵥ
    ρᴱᵍᵥ = epi_params.ρᴱᵍᵥ
    ρᴬᵍᵥ = epi_params.ρᴬᵍᵥ
    ρᴵᵍᵥ = epi_params.ρᴵᵍᵥ
    ρᴾᴴᵍᵥ = epi_params.ρᴾᴴᵍᵥ
    ρᴾᴰᵍᵥ = epi_params.ρᴾᴰᵍᵥ
    ρᴴᴿᵍᵥ = epi_params.ρᴴᴿᵍᵥ
    ρᴴᴰᵍᵥ = epi_params.ρᴴᴰᵍᵥ
    ρᴰᵍᵥ = epi_params.ρᴰᵍᵥ
    ρᴿᵍᵥ = epi_params.ρᴿᵍᵥ
    
    CHᵢᵍᵥ = epi_params.CHᵢᵍᵥ
    G = population.G
    M = population.M
    Nᵍ = population.Nᵍ
    nᵢᵍ = population.nᵢᵍ
    V = epi_params.V
    pᵍ_eff = population.pᵍ_eff
    C = population.C
    edgelist = population.edgelist
    Rᵢⱼ = population.Rᵢⱼ
    kᵍ_h = population.kᵍ_h
    kᵍ_hw = population.kᵍ_h .+ population.kᵍ_w
    
    # Total population
    N = sum(nᵢᵍ)

    # Intervention at time tᶜ
    if tᶜ == t
        pᵍ_eff[:] .= (1 - κ₀) .* population.pᵍ

        
        if (κ₀ != 0. )
            population.kᵍ_eff .= kᵍ_h * κ₀ .+ kᵍ_hw * (1 - δ) * (1 - κ₀)
            # elder keep home contacts during confinement
            population.kᵍ_eff[G] = kᵍ_h[G]
        end

        update_population_params!(population)
    end

    # Get P and compute Q
    compute_P!(Pᵢᵍᵥ, Sᵢᵍᵥ, pᵍ_eff, ρˢᵍᵥ, ρᴬᵍᵥ, ρᴵᵍᵥ,
               epi_params.Qᵢᵍ, population.nᵢᵍ_eff, population.mobilityᵍ,
               population.normᵍ, epi_params.βᴬ[1], epi_params.βᴵ[1],
               edgelist, Rᵢⱼ, C, M, G, V, t, rᵥ, kᵥ)
                
    
    # Compute τᵢᵍᵥ
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        
        for v in 1:V
            @simd for g in 1:G
                τᵢᵍᵥ[g, i, v] += Rᵢⱼ[indx_e] * Pᵢᵍᵥ[g, j, v]
            end
        end
    end
    
    
    # Compute vaccine priority distribution among ages and patches
    
    ϵᵢᵍ = optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                           ρˢᵍᵥ::Array{Float64, 4},
                                           nᵢᵍ::Array{Float64, 2},
                                           t::Int64)
            
    
    # Newly vaccinated people as a fraction of the subpopulation
    new_vaccinated = zeros(Float64, G, M)
    new_vaccinated .= ϵᵢᵍ ./ nᵢᵍ
    new_vaccinated[isnan.(new_vaccinated)] .= 0.

    # Update probabilities
    @inbounds for i in 1:M

        # Compute secure households
        CHᵢ = 0.0
        if tᶜ == t
            for g in 1:G
                aux = 0.0
                @simd for v in 1:V
                    aux += ρᴱᵍᵥ[g, i, t, v] + ρᴬᵍᵥ[g, i, t, v] + ρᴵᵍᵥ[g, i, t, v] 
                end
                CHᵢ += ( 1 - aux ) * population.nᵢᵍ[g, i]
            end
            CHᵢ = (1 - ϕ) * κ₀ * (CHᵢ / population.nᵢ[i]) ^ population.σ
        end
      
        
        # Update compartmental probabilities
        for g in 1:G
            if tᶜ == t
                ρˢᵍᵥ[g, i, t, :] .= ρˢᵍᵥ[g, i, t, :] .+ CHᵢᵍᵥ[g, i, t, :]
            end  
                  
            @simd for v in 1:V 
                
                # Infection probability
                Πᵢᵍᵥ = (1 - pᵍ_eff[g]) * Pᵢᵍᵥ[g, i, v] + pᵍ_eff[g] * τᵢᵍᵥ[g, i, v]

                # Pier: this is an ugly fix to avoid the problem of getting more vaccines that susceptibles
                # TO DO: incorporate this condition in the function optimal_vaccination_distribution
                if (v == 1) & ( (Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] + new_vaccinated[g, i]) > ρˢᵍᵥ[g, i, t, v] )
                    # print([g, i, t, v]  , "\n" )
                    new_vaccinated[g, i] = ρˢᵍᵥ[g, i, t, v] - Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] - 1e-10
                end

                #Epidemic compartments, where all states of vaccination are present
                ρˢᵍᵥ[g, i, t + 1, v] = ( 1 - Πᵢᵍᵥ ) * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] +
                    new_vaccinated[g, i] * ( [-1, +1, 0][v] ) +
                    Λ * ( [0, -1, +1][v] ) * ρˢᵍᵥ[g, i, t, 2] +  # The term inside the parentheses works as a if-then clause
                    Γ * ( [0 , 0, +1][v] ) * (ρᴿᵍᵥ[g, i, t, 1] + ρᴿᵍᵥ[g, i, t, 2] + ρᴿᵍᵥ[g, i, t, 3])
                
                ρᴱᵍᵥ[g, i, t + 1, v] = (1 - ηᵍ[g]) * ρᴱᵍᵥ[g, i, t, v] +
                    Πᵢᵍᵥ * (1 - CHᵢ) * ρˢᵍᵥ[g, i, t, v] 
        
                ρᴬᵍᵥ[g, i, t + 1, v] = (1 - αᵍ[g]) * ρᴬᵍᵥ[g, i, t, v] +
                    ηᵍ[g] * ρᴱᵍᵥ[g, i, t, v]

                ρᴵᵍᵥ[g, i, t + 1, v] = (1 - μᵍ[g]) * ρᴵᵍᵥ[g, i, t, v] +
                    αᵍ[g] * ρᴬᵍᵥ[g, i, t, v]

                ρᴾᴴᵍᵥ[g, i, t + 1, v] = (1 - λᵍ[g]) * ρᴾᴴᵍᵥ[g, i, t, v] +
                    μᵍ[g] * (1 - θᵍᵥ[g, v]) * γᵍᵥ[g, v] *  ρᴵᵍᵥ[g, i, t, v]

                ρᴾᴰᵍᵥ[g, i, t + 1, v] = (1 - ζᵍ[g]) * ρᴾᴰᵍᵥ[g, i, t, v] +
                    μᵍ[g] * θᵍᵥ[g, v] * ρᴵᵍᵥ[g, i, t, v]

                ρᴴᴿᵍᵥ[g, i, t + 1, v] = (1 - χᵍ[g]) * ρᴴᴿᵍᵥ[g, i, t, v] +
                    λᵍ[g] * (1 - ωᵍᵥ[g, v] ) * ρᴾᴴᵍᵥ[g, i, t, v]

                ρᴴᴰᵍᵥ[g, i, t + 1, v] = (1 - ψᵍ[g]) * ρᴴᴰᵍᵥ[g, i, t, v] +
                    λᵍ[g] * ωᵍᵥ[g, v] * ρᴾᴴᵍᵥ[g, i, t, v]

                ρᴿᵍᵥ[g, i, t + 1, v] = ρᴿᵍᵥ[g, i, t, v] + χᵍ[g] * ρᴴᴿᵍᵥ[g, i, t, v] +
                    μᵍ[g] * (1 - θᵍᵥ[g, v]) * (1 - γᵍᵥ[g, v]) * ρᴵᵍᵥ[g, i , t, v] -
                    Γ * ρᴿᵍᵥ[g, i, t, v] 

                ρᴰᵍᵥ[g, i, t + 1, v] = ρᴰᵍᵥ[g, i, t, v] + ζᵍ[g] * ρᴾᴰᵍᵥ[g, i, t, v] +
                    ψᵍ[g] * ρᴴᴰᵍᵥ[g, i, t, v]
                
            end

            if tᶜ == t
                for v in 1:V
                    aux = ρˢᵍᵥ[g, i, t, v]
                    ρˢᵍᵥ[g, i, t, v] -= CHᵢᵍᵥ[g, i, t, v] 
                    CHᵢᵍᵥ[g, i, t + 1 , v] = CHᵢ * aux
                end
            end  
            
            # Reset values
            τᵢᵍᵥ[g, i, :] .= 0.
	        # this should be one, based on the intial value provided in run_epidemic_spreading_mmca
            Pᵢᵍᵥ[g, i, :] .= 0. 
        end
    end

end


"""
    compute_P!(Pᵢᵍᵥ::Array{Float64, 3},
               Sᵢᵍᵥ::Array{Float64, 3},
               pᵍ_eff::Array{Float64, 1},
               ρˢᵍᵥ::Array{Float64, 4},
               ρᴬᵍᵥ::Array{Float64, 4},
               ρᴵᵍᵥ::Array{Float64, 4},
               Qᵢᵍ::Array{Float64, 3},
               nᵢᵍ_eff::Array{Float64, 2},
               mobilityᵍ::Array{Float64, 2},
               normᵍ::Array{Float64, 2},
               βᴬ::Float64,
               βᴵ::Float64,
               edgelist::Array{Int64, 2},
               Rᵢⱼ::Array{Float64, 1},
               C::Array{Float64, 2},
               M::Int64,
               G::Int64,
               V::Int64,
               t::Int64,
               rᵥ::Array{Float64, 1},
               kᵥ::Array{Float64, 1})

Compute ``P_i^g_v(t)`` and ``Q_i^g(t)`` as described in the referenced paper. The first quantity is needed to compute the infection probability and the second for the basig reproduction number
"""

function compute_P!(Pᵢᵍᵥ::Array{Float64, 3},
                    Sᵢᵍᵥ::Array{Float64, 3},
                    pᵍ_eff::Array{Float64, 1},
                    ρˢᵍᵥ::Array{Float64, 4},
                    ρᴬᵍᵥ::Array{Float64, 4},
                    ρᴵᵍᵥ::Array{Float64, 4},
                    Qᵢᵍ::Array{Float64, 3},
                    nᵢᵍ_eff::Array{Float64, 2},
                    mobilityᵍ::Array{Float64, 2},
                    normᵍ::Array{Float64, 2},
                    βᴬ::Float64,
                    βᴵ::Float64,
                    edgelist::Array{Int64, 2},
                    Rᵢⱼ::Array{Float64, 1},
                    C::Array{Float64, 2},
                    M::Int64,
                    G::Int64,
                    V::Int64,
                    t::Int64,
                    rᵥ::Array{Float64, 1},
                    kᵥ::Array{Float64, 1})
    

    # Init. aux variables
    Sᵢᵍᵥ = zeros(G, M, V)
    Pᵢᴬᵍᵥ = zeros(G, M, V)
    Pᵢᴵᵍᵥ = zeros(G, M, V)
    nˢᵍᵥ_ij = zeros(V)
    
    nᴬᵍᵥ_ij = zeros(G, M, V)
    nᴵᵍᵥ_ij = zeros(G, M, V)

    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2] # i->j

        # Get effective S, A and I
        for g in 1:G
            nˢᵍᵥ_ij[:] .= ( ρˢᵍᵥ[g, i, t, :] .* (1 .- rᵥ[:]) ) * mobilityᵍ[g, indx_e]
            Sᵢᵍᵥ[g, j, :] .=  Sᵢᵍᵥ[g, j, :] .+ nˢᵍᵥ_ij[:] / nᵢᵍ_eff[g, j]
            @simd for h in 1:G
                nᴬᵍᵥ_ij[g, j, :] .= nᴬᵍᵥ_ij[g, j, :] .+ ρᴬᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
                nᴵᵍᵥ_ij[g, j, :] .= nᴵᵍᵥ_ij[g, j, :] .+ ρᴵᵍᵥ[h, i, t, :] .* (C[g, h] * mobilityᵍ[h, indx_e] ./ nᵢᵍ_eff[h, j] .* ones(V) )
            end
        end
    end


    # @inbounds for indx_e in 1:length(Rᵢⱼ)
    #     i = edgelist[indx_e, 1]
    #     j = edgelist[indx_e, 2]

    #     # Get effective S, A and I
    #     for g in 1:G
    #         nˢᵍ_ij = ρˢᵍ[g, i, t] * mobilityᵍ[g, indx_e]
    #         Sᵢᵍ[g, j] +=  nˢᵍ_ij / nᵢᵍ_eff[g, j]
    #         @simd for h in 1:G
    #             nᴬᵍ_ij = ρᴬᵍ[h, i, t] * mobilityᵍ[h, indx_e]
    #             nᴵᵍ_ij = ρᴵᵍ[h, i, t] * mobilityᵍ[h, indx_e]
    #             Pᵢᴬᵍ[g, j] += C[g, h] * nᴬᵍ_ij / nᵢᵍ_eff[h, j]
    #             Pᵢᴵᵍ[g, j] += C[g, h] * nᴵᵍ_ij / nᵢᵍ_eff[h, j]
    #         end
    #     end
    # end

    # Get P and effective ρ
    @inbounds for i in 1:M
        for v in 1:V
            for g in 1:G
                aux = 1.0
                @simd for w in 1:V
                    aux = aux * (1 - βᴬ*(1 - rᵥ[v])*(1 - kᵥ[w]) )^(normᵍ[g, i] * nᴬᵍᵥ_ij[g, i, w]) *
                                (1 - βᴵ*(1 - rᵥ[v])*(1 - kᵥ[w]) )^(normᵍ[g, i] * nᴵᵍᵥ_ij[g, i, w]) 
                end
                Pᵢᵍᵥ[g, i, v] = 1 - aux
            end
        end
    end
    
    # Compute Q to get the effective R
    @inbounds for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2] # i->j
        
        for g in 1:G     
            for v in 1:V
                @simd for h in 1:G
                    Qᵢᵍ[g, i, t] += normᵍ[g, j] * C[g, h] * Sᵢᵍᵥ[h, j, v] *
                    (pᵍ_eff[g] * Rᵢⱼ[indx_e] + (1 - pᵍ_eff[g]) * (i == j ? 1. : 0.))
                end
            end
        end
    end

end



"""
    update_population_params!(population::Population_Params)

Update population parameters, computing the effective populations and the
normalization parameter z if p and k are modified.

# Arguments
  - `population::Population_Params`: Structure that contains all the parameters related with the population.
"""
function update_population_params!(population::Population_Params)

    # Reset effective population
    population.nᵢᵍ_eff[:,:]    .= (1 .- population.pᵍ_eff) .* population.nᵢᵍ
    population.nᵢ_eff[:]       .= zeros(Float64, population.M)
    population.mobilityᵍ[:, :] .= 0.

    # Init. normalization vector
    population.zᵍ .= 0.
    population.normᵍ[:, :] .= 0.

    # Compute effective population
    compute_effective_population!(population.nᵢᵍ_eff,
                                  population.nᵢ_eff,
                                  population.nᵢᵍ,
                                  population.Nᵍ,
                                  population.mobilityᵍ,
                                  population.kᵍ_eff,
                                  population.zᵍ,
                                  population.normᵍ,
                                  population.ξ,
                                  population.pᵍ_eff,
                                  population.sᵢ,
                                  population.edgelist,
                                  population.Rᵢⱼ,
                                  population.M,
                                  population.G)
end


"""
    compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                  nᵢ_eff::Array{Float64, 1},
                                  nᵢᵍ::Array{Float64, 2},
                                  Nᵍ::Array{Int64, 1},
                                  mobilityᵍ::Array{Float64, 2},
                                  kᵍ_eff::Array{Float64, 1},
                                  zᵍ::Array{Float64, 1},
                                  normᵍ::Array{Float64, 2},
                                  ξ::Float64,
                                  pᵍ_eff::Array{Float64, 1},
                                  sᵢ::Array{Float64, 1},
                                  edgelist::Array{Int64, 2},
                                  Rᵢⱼ::Array{Float64, 1},
                                  M::Int64,
                                  G::Int64)

Compute the effective population at each patch.
"""
function compute_effective_population!(nᵢᵍ_eff::Array{Float64, 2},
                                       nᵢ_eff::Array{Float64, 1},
                                       nᵢᵍ::Array{Float64, 2},
                                       Nᵍ::Array{Int64, 1},
                                       mobilityᵍ::Array{Float64, 2},
                                       kᵍ_eff::Array{Float64, 1},
                                       zᵍ::Array{Float64, 1},
                                       normᵍ::Array{Float64, 2},
                                       ξ::Float64,
                                       pᵍ_eff::Array{Float64, 1},
                                       sᵢ::Array{Float64, 1},
                                       edgelist::Array{Int64, 2},
                                       Rᵢⱼ::Array{Float64, 1},
                                       M::Int64,
                                       G::Int64)

    # Compute the effective population for each strata and patch
    for indx_e in 1:size(edgelist)[1]
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]

        for g in 1:G
            nᵢᵍ_eff[g, j] += pᵍ_eff[g] * Rᵢⱼ[indx_e] * nᵢᵍ[g, i]
        end
    end

    # Compute the aggregated effective population
    for i in 1:M
        for g in 1:G
            nᵢ_eff[i] += nᵢᵍ_eff[g, i]
        end
    end

    # Correction for populations without effective population
    for i in 1:M
        if nᵢ_eff[i] == 0.
            nᵢ_eff[i] = 1e-7
        end
        for g in 1:G
            if nᵢᵍ_eff[g, i] == 0.
                nᵢᵍ_eff[g, i] = 1e-7
            end
        end
    end

    # Compute the normalization vector
    for g in 1:G
        zᵍ[g] = kᵍ_eff[g] * Nᵍ[g] /
            sum((2 .- exp.(-ξ .* nᵢ_eff ./ sᵢ)) .* nᵢᵍ_eff[g, :]);
    end

    # Update the precomuted matrices
    for indx_e in 1:length(Rᵢⱼ)
        i = edgelist[indx_e, 1]
        j = edgelist[indx_e, 2]
        for g in 1:G
            mobilityᵍ[g, indx_e] = nᵢᵍ[g, i] *
                ((1 - pᵍ_eff[g]) * (i == j ? 1. : 0.) + pᵍ_eff[g] * Rᵢⱼ[indx_e] )
        end
    end

    for i in 1:M
        for g in 1:G
            normᵍ[g, i] = zᵍ[g] * (2 - exp(-ξ * nᵢ_eff[i] / sᵢ[i]))
        end
    end
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params;
                                 tᶜ::Int64 = -1,
                                 tᵛ::Int64 = -1,
                                 κ₀::Float64 = 0.0,
                                 ϕ::Float64 = 1.0,
                                 δ::Float64 = 0.0,
                                 ϵᵍ::Array{Float64, 1} = [0., 0., 0.],
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It also provides, through optional arguments,
the application of a containmnet or a vaccination campaign on a specific date.

# Arguments
    - `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
    and the epidemic spreading information.
    - `population::Population_Params`: Structure that contains all the parameters
    related with the population.

## Optional
    - `tᶜ::Int64 = -1`: Timestep of application of containment, or out of timesteps range
    value for no containment.
    - `tᵛ::Int64 = -1`: Timestep of application of vaccination.
    - `κ⁰::Float64 = 0.0`: Mobility reduction.
    - `ϕ::Float64 = 1.0`: Permeability of confined households.
    - `δ::Float64 = 0.0`: Social Distancing.
    - `ϵᵍ::Array{Float64, 1} = [0., 0., 0.]`: Number of vaccines for each age group
    - `t₀::Int64 = 1`: Initial timestep.
    - `verbose::Bool = false`: If `true`, prints useful information about the
    evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      npi_params::NPI_Params;
                                      tᵛ::Int64 = -1,
                                      ϵᵍ::Array{Float64, 1} = [0., 0., 0.],
                                      t₀::Int64 = 1,
                                      verbose::Bool = false)
    G = population.G
    M = population.M
    V = epi_params.V

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍᵥ = zeros(Float64, G, M, V)
    Pᵢᵍᵥ = zeros(Float64, G, M, V)

    # Auxiliar arrays to compute P (avoid the allocation of additional memory)
    Sᵢᵍᵥ = zeros(Float64, G, M, V)
    
    
    run_epidemic_spreading_mmca!(epi_params, population, npi_params, [tᵛ], reshape(ϵᵍ, (3,1)) , t₀ = t₀, verbose = verbose)
end


"""
    run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                 population::Population_Params,
                                 tᶜs::Array{Int64, 1},
                                 tᵛs::Array{Int64, 1},
                                 κ₀s::Array{Float64, 1},
                                 ϕs::Array{Float64, 1},
                                 δs::Array{Float64, 1},
                                 ϵᵍs::Array{Float64, 2};
                                 t₀::Int64 = 1,
                                 verbose::Bool = false)

Computes the evolution of the epidemic spreading over time, updating the
variables stored in epi_params. It provides the option of the application
of multiple different containmnets or vaccination campaigns at specific dates.

# Arguments
    - `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
    and the epidemic spreading information.
    - `population::Population_Params`: Structure that contains all the parameters
    related with the population.
    - `tᶜs::Array{Int64, 1}`: List of timesteps of application of containments.
    - `tᵛs::Int64 = -1`: Timestep of application of vaccination.
    - `κ⁰s::Array{Float64, 1}`: List of mobility reductions.
    - `ϕs::Array{Float64, 1}`: List of permeabilities of confined households.
    - `δs::Array{Float64, 1}`: List of social distancings.
    - `ϵᵍs::Array{Float64, 2}`: List of dosis per age group of each time period

## Optional
    - `t₀::Int64 = 1`: Initial timestep.
    - `verbose::Bool = false`: If `true`, prints useful information about the evolution of the epidemic process.
"""
function run_epidemic_spreading_mmca!(epi_params::Epidemic_Params,
                                      population::Population_Params,
                                      npi_params::NPI_Params,
                                      tᵛs::Array{Int64, 1},
                                      ϵᵍs::Array{Float64, 2};
                                      t₀::Int64 = 1,
                                      verbose::Bool = false)

    G = population.G
    M = population.M
    T = epi_params.T
    V = epi_params.V

    # Initialize τᵢ (Π = (1 - p) P + pτ) and Pᵢ for markov chain
    τᵢᵍᵥ = zeros(Float64, G, M, V)
    Pᵢᵍᵥ = ones(Float64, G, M, V)

    Sᵢᵍᵥ = zeros(Float64, G, M, V)

    κ₀s = npi_params.κ₀s
     ϕs = npi_params.ϕs
    δs = npi_params.δs
    tᶜs = npi_params.tᶜs

    # Initial state
    if verbose
        print_status(epi_params, population, t₀)
    end

    i = 1 # counter for containment
    j = 1 # counter for vaccination

    ## Start loop for time evoluiton
    @inbounds for t in t₀:(T - 1)
        
        update_prob!(Pᵢᵍᵥ, Sᵢᵍᵥ, τᵢᵍᵥ, epi_params, population,
                        κ₀s[i], ϕs[i], δs[i], ϵᵍs[:, j], t, tᶜs[i], tᵛs[j])
        
        if t == tᶜs[i] && i < length(tᶜs)
            i += 1
        end

        if t == tᵛs[j] && j < length(tᵛs)
            j += 1
        end
        
        # To avoid negative compartments
        only_positive = true
        for v in 1:V
            only_positive = only_positive & all(epi_params.ρˢᵍᵥ[:, :, t, v] .>= 0.0) & all(epi_params.ρˢᵍᵥ[:, :, t, v] .<= 1.0)
        end
        
        if !only_positive
            @printf("ATTENZIONE: I suscettibili sono meno di 0 o più di 1\n")
            return
        end

        if verbose
            print_status(epi_params, population, t + 1)
        end
    end
end
