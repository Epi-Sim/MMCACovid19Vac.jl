### ----------------------------------------------------------------------------
### COMPUTE R FUNCTIONS
### ----------------------------------------------------------------------------

"""
    compute_R_eff(epi_params::Epidemic_Params,
                  population::Population_Params;
                  τ::Int64 = 21)

Compute the effective reproduction number R.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.

# Optional

- `τ::Int64 = 21`: kernel length.

# Return

Tuple formed by:
- DataFrame containing information about the evolution of the effective
  reproduction number R for each strata and patch.
- DataFrame containing information about the evolution of the total effective
  reproduction number R.
"""
function compute_R_eff(epi_params::Epidemic_Params,
                       population::Population_Params,
                       τ::Int64 = 21)

    M = population.M
    G = population.G
    T = epi_params.T - 1 # because we have to cut out the case t=1
    V = epi_params.V
    
    βᴬ = epi_params.βᴬ[1]
    βᴵ = epi_params.βᴵ[1]

    ηᵍ = epi_params.ηᵍ
    αᵍ = epi_params.αᵍ
    μᵍ = epi_params.μᵍ
    
    rᵥ = epi_params.rᵥ
    kᵥ = epi_params.kᵥ
    
    ρˢᵍᵥ = epi_params.ρˢᵍᵥ

    # Setup the kernel
    tʷ = 0:(τ - 1)
    wᵍᵥ = ones(Float64, G, τ, V)

    # Initialize results
    Rᵢᵍ_eff = DataFrame()
    Rᵢᵍ_eff.strata = repeat(1:G, outer = (T - τ) * M * V)
    Rᵢᵍ_eff.vaccine = repeat(1:V, inner = G, outer = (T - τ) * M)
    Rᵢᵍ_eff.patch = repeat(1:M, inner = G * V, outer = (T - τ))
    Rᵢᵍ_eff.time = repeat(1:(T - τ), inner = G * M * V)

    R_eff = DataFrame()
    R_eff.time = 1:(T - τ)
    Δρ = ones(G, M, T - τ, V)

    # Compute R
    Rᵢᵍᵥ = zeros(Float64, G, M, (T - τ), V)
    Rᵢᵍ = ones(Float64, G, M, (T - τ) )
    Rᵍ = ones(Float64, G, (T - τ) )
    Rᵍᵥ= ones(Float64, G, (T - τ), V )
    Rᵥ = ones(Float64, (T - τ), V )
    R = zeros(T - τ)
    for t in 1:(T - τ)
        for v in 1:V
    # This calculation would hold only without reinfection and waning immunity
    #   @simd for s in tʷ
    #       @. Rᵢᵍᵥ[:, :, t, v] += epi_params.Qᵢᵍ[:, :, t + s + 2] * wᵍᵥ[:, s + 1, v] 
    #       # the +2 instead of +1 is because time t corresponds to the time t+1 in the epi_params
    #   end
            @simd for g in 1:G
                @. Rᵢᵍᵥ[g, :, t, v] = epi_params.Qᵢᵍ[g, :, t] * (1 - kᵥ[v]) * (βᴬ/αᵍ[g] + βᴵ/μᵍ[g] )
            end
        end

        # Notice that Δρ of time t corresponds to the time t+1 in the epi_params
        # @. Δρ[:, :, t, :] = ρˢᵍᵥ[:, :, t, :] - ρˢᵍᵥ[:, :, t + 1, :]
        # R[t] = sum(Rᵢᵍᵥ[:, :, t, :] .* Δρ[:, :, t, :] .* population.nᵢᵍ[:, :]) / sum( Δρ[:, :, t, :] .* population.nᵢᵍ[:, :] )
        
        Rᵍᵥ[:, t, :] = sum( Rᵢᵍᵥ[:, :, t, :] , dims = 2 )[:, 1, :] / M
        Rᵥ[t, :] = sum( Rᵍᵥ[:, t, :] , dims = 1 )[1, :] / G 
        Rᵢᵍ[:, :, t] = sum( Rᵢᵍᵥ[:, :, t, :] , dims = 4 )[:,:,1] / V
        Rᵍ[:, t] = sum( Rᵢᵍ[:, :, t] , dims = 2 )[:,1] / M 
        R[t] = sum(Rᵢᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :]) / population.N
        
    end

    Rᵢᵍ_eff.R_eff = reshape(Rᵢᵍᵥ, G * M * (T - τ) * V)
    R_eff.R_eff = R

    return Rᵢᵍ_eff, Rᵍ, R_eff, Rᵥ, Rᵍᵥ
end


"""
    Function whose purpose is to find all the local maxima in a vector. The output is a list of two vectors, one containing all the heights 
    of said maxima, the other containing all the positions
"""

function maxima(v)
    h = []
    p = []
    idx = 1
    for i in 2:(length(v)-2)
        condition = (v[i-1] < v[i]) & (v[i] > v[i+1])
        # condition2 = (v[maximum([i-50, 1])] < v[i]) | (v[i] > v[minimum([length(v), i+50])])
        if condition #& condition2
           h = append!( h, v[i] )
           p = append!( p, i)
           idx = idx + 1 
        end
    end
    return (height = h ,
            position = p)
end


"""
    optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                     ρˢᵍᵥ::Array{Float64, 4},
                                     nᵢᵍ::Array{Float64, 2},
                                     t::Int64)

Computes the number of vaccines that should be distributed among spatial patches and 
age strata taking into account where the majority of the unvaccinated susceptibles 
are and what age stratum should have the priority.

# Arguments

- ϵᵍ : the absolute number of vaccines at our disposal per day, divided among age strata
- ρˢᵍᵥ: fraction of susceptible individuals
- nᵢᵍ: absolute number of people per age in each age strata (rows) and patches (colums)
- t: current day

"""

function optimal_vaccination_distribution(ϵᵍ::Array{Float64, 1},
                                          ρˢᵍᵥ::Array{Float64, 4},
                                          nᵢᵍ::Array{Float64, 2},
                                          t::Int64)
    Nᵥ = sum(ϵᵍ) # Total number of vaccines
    (G, M, T, V) = size(ρˢᵍᵥ)
    (G, M) = size(nᵢᵍ)

    if Nᵥ == 0
        @debug "No vaccination"
        return zeros(G, M)
    end
    
    
    if any(ϵᵍ .< 0)
        @error "COMPUTING ERROR: Number of dosis is negative"
        return
    end

    only_positive = true
    for v in 1:V
        only_positive = only_positive & all(ρˢᵍᵥ[:, :, t, v] .>= 0.0) & all(ρˢᵍᵥ[:, :, t, v] .<= 1.0)
    end
    
    if !only_positive
        @error "COMPUTING ERROR: Fraction of susceptible is not between 0 and 1"
        return
    end
    
    ###############################
    
    if ( ( Nᵥ != 0) & (sum(ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ) > Nᵥ) ) 
         
        # Define a matrix that gives you the priority of each patch and age group...
        priority_ϵ =  nᵢᵍ .* ( reshape(repeat(ϵᵍ, M), (G,M) ) )
        priority_ϵ = priority_ϵ / (sum(priority_ϵ) == 0 ? 1 : sum(priority_ϵ) )
        # ... and use the priority matrix to define how many dosis each subgroup get
        ϵᵢᵍ = Nᵥ * priority_ϵ
        
        # Define index that tells you if and where there are more susceptibles than vaccines
        idx_ϵ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .- ϵᵢᵍ

        # Redistribution of vaccines: If there is one location and age group that has more vaccines than susceptibles the number 
        # of vaccines in that compartment is set equal to the number of susceptibles and the spare dosis are restributed among the others
        while ( !prod(idx_ϵ .>= 0 ) )   
            ϵᵢᵍ = ϵᵢᵍ .* (idx_ϵ .> 0)
            ϵᵢᵍ .* (idx_ϵ .<= 0) .= ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ .* (idx_ϵ .<= 0)
            Nᵥ_new = Nᵥ - sum(nᵢᵍ .* (idx_ϵ .<= 0) )
            
            # Redifine priority levels
            priority_ϵ =  nᵢᵍ .* ϵᵢᵍ 
            priority_ϵ = priority_ϵ / (sum(priority_ϵ) == 0 ? 1 : sum(priority_ϵ) )
            ϵᵢᵍ .* (idx_ϵ .> 0) .= Nᵥ_new * priority_ϵ
            # ϵᵢᵍ = ϵᵢᵍ / sum(ϵᵢᵍ) * ( Nᵥ - sum(nᵢᵍ[idx_ϵ .<= 0])  )
            
            idx_ϵ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ  .-  ϵᵢᵍ[:, :] 
        end
        
    elseif ( (Nᵥ != 0) & (sum(ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ) <= Nᵥ)  )
        ϵᵢᵍ = ρˢᵍᵥ[:,:,t,1] .* nᵢᵍ
    else
        ϵᵢᵍ = zeros(G, M)
    end
    
    return ϵᵢᵍ

end




"""
Get the time series of different indicators of the epidemic:
 - Infected = I + A
 - Cases = I + A + PD + PH + HD + HR
 - Icus = HR + HD
 - Deaths = D
 - Vaccinated = All the vaccinated compartments
 - Daily cases
"""

function time_series(epi_params::Epidemic_Params,
                      population::Population_Params)
    
    return (infected = sum((epi_params.ρᴵᵍᵥ[:, :, :, :] .+ 
                            epi_params.ρᴬᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            cases    = sum((epi_params.ρᴵᵍᵥ[:, :, :, :] .+ 
                            epi_params.ρᴬᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴾᴰᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴾᴴᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴰᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴿᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            icus     = sum((epi_params.ρᴴᴿᵍᵥ[:, :, :, :] .+
                            epi_params.ρᴴᴰᵍᵥ[:, :, :, :]) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            deaths   = sum(epi_params.ρᴰᵍᵥ[:, :, :, :] .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1],
        
            vaccinated = sum((epi_params.ρˢᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴾᴰᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴱᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴬᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴵᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴾᴴᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴴᴰᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴴᴿᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴿᵍᵥ[:, :, :, 2:3] .+
                              epi_params.ρᴰᵍᵥ[:, :, :, 2:3] ) .* population.nᵢᵍ[:, :], dims=(1,2,4) )[1,1,:,1]
        
        )
end


