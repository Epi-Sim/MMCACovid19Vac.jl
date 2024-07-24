function create_default_epi_params()
    epiparams_dict = Dict()
    epiparams_dict["scale_β"] = 0.51
    epiparams_dict["βᴬ"] = 0.046053
    epiparams_dict["βᴵ"] = 0.0903
    epiparams_dict["ηᵍ"] = [0.2747252747252747, 0.2747252747252747, 0.2747252747252747]
    epiparams_dict["αᵍ"] = [0.26595744680851063, 0.641025641025641, 0.641025641025641]
    epiparams_dict["μᵍ"] = [1.0, 0.3125, 0.3125]
    epiparams_dict["θᵍ"] = [0.0, 0.0, 0.0]
    epiparams_dict["γᵍ"] = [0.003, 0.01, 0.08]
    epiparams_dict["ζᵍ"] = [0.12820512820512822, 0.12820512820512822, 0.12820512820512822]
    epiparams_dict["λᵍ"] = [1.0, 1.0, 1.0]
    epiparams_dict["ωᵍ"] = [0.0, 0.04, 0.3]
    epiparams_dict["ψᵍ"] = [0.14285714285714285, 0.14285714285714285, 0.14285714285714285]
    epiparams_dict["χᵍ"] = [0.047619047619047616, 0.047619047619047616, 0.047619047619047616]
    epiparams_dict["Λ"] = 0.02
    epiparams_dict["Γ"] = 0.01
    epiparams_dict["rᵥ"] = [0.0, 0.6, 0.0]
    epiparams_dict["kᵥ"] = [0.0, 0.4, 0.0]
    epiparams_dict["risk_reduction_dd"] = 0.0
    epiparams_dict["risk_reduction_h"] = 0.1
    epiparams_dict["risk_reduction_d"] = 0.05
    
    return epiparams_dict
end

function create_default_population_params()
    population = Dict()
    populaiont["age_labels"] = ["Y", "M", "O"]
    populaiont["C"] = [ 0.598  0.38486 0.01714 ;
                        0.244  0.721   0.0353;
                        0.1919 0.5705  0.2376
                      ]
    population["kᵍ"] = [11.8, 13.3, 6.76]
    population["kᵍ_h"] = [3.15, 3.17, 3.28]
    population["kᵍ_w"] = [1.72, 5.18, 0.0]
    population["pᵍ"] = [0.0, 1.0, 0.00]
    population["ξ"] = 0.01
    population["σ"] = 2.5
    
    return population
end

function create_default_vacparameters()
    vacparams_dict = Dict()
    vacparams_dict["ϵᵍ"] = [0.1 , 0.4 , 0.5]
    vacparams_dict["percentage_of_vacc_per_day"] = 0.005
    vacparams_dict["start_vacc"] = 2
    vacparams_dict["dur_vacc"] = 8
    vacparams_dict["are_there_vaccines"] = false

    return vacparams_dict
end

function create_default_npiparameters()
    npiparams_dict = Dict()
    ## It's important that the default parameters are those of the absence of 
    ## lockdowns, because they are the one the code refers to if the key "are_there_npi" = false
    npiparams_dict["κ₀s"] = [0.0]
    npiparams_dict["ϕs"] = [1.0]
    npiparams_dict["δs"] = [0.0]
    npiparams_dict["tᶜs"] =  [1]
    npiparams_dict["are_there_npi"] = true

    return npiparams_dict
end

function update_config!(config, cmd_line_args)
    # Define dictionary containing epidemic parameters
    if !haskey(config, "epidemic_params")
        config["epidemic_params"] = create_default_epi_params()
    end

    if !haskey(config, "population_params")
        config["population_params"] = create_default_population_params()
    end

    # Define dictionary containing vaccination parameters
    if !haskey(config, "vaccination")
        config["vaccination"] = create_default_vacparameters()
    end

    # Define dictionary containing npi parameters
    if (!haskey(config, "NPI") | !config["NPI"]["are_there_npi"] )
        config["NPI"] = create_default_npiparameters()
    end

    # overwrite config with command line
    if cmd_line_args["start-date"] !== nothing
        config["simulation"]["first_day_simulation"] = cmd_line_args["start-date"]
    end
    if cmd_line_args["end-date"] !== nothing
        config["simulation"]["last_day_simulation"] = cmd_line_args["end-date"]
    end
    if cmd_line_args["export-compartments-time-t"] !== nothing
        config["simulation"]["export_compartments_time_t"] = cmd_line_args["export-compartments-time-t"]
    end
    if cmd_line_args["export-compartments-full"] == true
        config["simulation"]["export_compartments_full"] = true
    end

    if cmd_line_args["initial-conditions"] !== nothing
        config["data"]["initial-conditions"] = cmd_line_args["initial-conditions"]
    end

    nothing
end

"""
save_simulation_hdf5(epi_params::Epidemic_Params,
                         population::Population_Params,
                         output_fname::String;
                         export_time_t = -1)

    Save the full simulations.

    # Arguments

    - `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
    and the epidemic spreading information.
    - `population::Population_Params`: Structure that contains all the parameters
    related with the population.
    - `output_fname::String`: Output filename.

    ## Optional

    - `export_time_t = -1`: Time step to ve saved instead of the full simulation.
"""
function save_simulation_hdf5(epi_params::Epidemic_Params, 
                              population::Population_Params,
                              output_fname;
                              export_time_t = -1)

    G = population.G
    M = population.M
    T = epi_params.T
    V = epi_params.V
    N = epi_params.NumComps

    compartments = zeros(Float64, G, M, T, V, N);
    compartments[:, :, :, :, 1]  .= epi_params.ρˢᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 2]  .= epi_params.ρᴱᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 3]  .= epi_params.ρᴬᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 4]  .= epi_params.ρᴵᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 5]  .= epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 6]  .= epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 7]  .= epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 8]  .= epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 9]  .= epi_params.ρᴿᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 10] .= epi_params.ρᴰᵍᵥ .* population.nᵢᵍ
    if export_time_t > 0
        h5open(output_fname, "w") do file
            write(file, "data", compartments[:,:,export_time_t,:,:])
        end
    else
        h5open(output_fname, "w") do file
            write(file, "data", compartments[:,:,:,:,:])
        end
    end
end


"""
save_simulation_netCDF(epi_params::Epidemic_Params, 
                                    population::Population_Params,
                                    output_fname::String;
                                    G_coords= nothing,
                                    M_coords = nothing,
                                    T_coords = nothing
                                    )

    Save the full simulations.

    # Arguments

    - `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
    and the epidemic spreading information.
    - `population::Population_Params`: Structure that contains all the parameters
    related with the population.
    - `output_fname::String`: Output filename.

    ## Optional
    - `G_coords = nothing`: Array::{String} of size G containing the labels for age strata
    - `M_coords = nothing`: Array::{String} of size M containing the labels for the patches
    - `T_coords = nothing`: Array::{String} of size t containing the labels for the time (dates)
    - `export_time_t = -1`: Time step to ve saved instead of the full simulation.
"""
function save_simulation_netCDF( epi_params::Epidemic_Params, 
                                 population::Population_Params,
                                 output_fname::String;
                                 G_coords = nothing,
                                 M_coords = nothing,
                                 T_coords = nothing,
                                 V_coords = nothing
                                )
    G = population.G
    M = population.M
    T = epi_params.T
    V = epi_params.V
    S = epi_params.NumComps
    S_coords = epi_params.CompLabels
    V_coords = epi_params.VaccLabels

    if isnothing(G_coords)
        G_coords = collect(1:G)
    end
    if isnothing(M_coords)
        M_coords = collect(1:M)
    end
    if isnothing(T_coords)
        T_coords = collect(1:T) 
    end

    compartments = zeros(Float64, G, M, T, V, S);

    # Adding 
    compartments[:, :, :, :, 1]  .= (epi_params.ρˢᵍᵥ + epi_params.CHᵢᵍᵥ) .* population.nᵢᵍ
    compartments[:, :, :, :, 2]  .= epi_params.ρᴱᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 3]  .= epi_params.ρᴬᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 4]  .= epi_params.ρᴵᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 5]  .= epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 6]  .= epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 7]  .= epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 8]  .= epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 9]  .= epi_params.ρᴿᵍᵥ .* population.nᵢᵍ
    compartments[:, :, :, :, 10] .= epi_params.ρᴰᵍᵥ .* population.nᵢᵍ
    isfile(output_fname) && rm(output_fname)

    nccreate(output_fname, "data", "G", G_coords, "M", M_coords, "T", T_coords, "V", V_coords, "epi_states", S_coords)
    ncwrite(compartments, output_fname, "data")

end


"""
    store_R_eff(epi_params::Epidemic_Params,
                population::Population_Params,
                suffix::String,
                folder::String;
                τ::Int64 = 21)

Compute and store the effective reproduction number R.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path where the results will be stored.

# Optional

- `τ::Int64 = 21`: kernel length.
"""
function store_R_eff(epi_params::Epidemic_Params,
                     population::Population_Params,
                     suffix::String,
                     folder::String;
                     τ::Int64 = 21)

    Rᵢᵍ_eff, R_eff = compute_R_eff(epi_params, population, τ)

    # Write the results
    CSV.write(@sprintf("%s/output_Reff_%s.csv", folder, suffix), Rᵢᵍ_eff)
    CSV.write(@sprintf("%s/output_Reff_total_%s.csv", folder, suffix), R_eff)
end


### ----------------------------------------------------------------------------
### OUTPUT FUNCTIONS
### ----------------------------------------------------------------------------

"""
    store_compartment(epi_params::Epidemic_Params,
                      population::Population_Params,
                      compartment::Char,
                      suffix::String,
                      folder::String)

Store the evolution of the given epidemic compartment for each strata and patch.

# Arguments

- `epi_params::Epidemic_Params`: Structure that contains all epidemic parameters
  and the epidemic spreading information.
- `population::Population_Params`: Structure that contains all the parameters
  related with the population.
- `compartment::String`: String indicating the compartment, one of: `"S"`,
  `"E"`, `"A"`, `"I"`, `"PH"`, `"PD"`, `"HR"`, `"HD"`, `"D"`, `"R"`.
- `suffix::String`: String used to identify the experiment.
- `folder::String`: String containing the path to the folder where the results
  will be stored.
"""
function store_compartment(epi_params::Epidemic_Params,
                           population::Population_Params,
                           compartment::String,
                           suffix::String,
                           folder::String)

    M = population.M
    G = population.G
    T = epi_params.T
    V = epi_params.V

    # Init. dataframe
    df = DataFrame()
    df.strata = repeat(1:G, outer = T * M)
    df.patch = repeat(1:M, inner = G, outer = T)
    df.time = repeat(1:T, inner = G * M)

    # Store number of cases
    if compartment == "S"
        df.cases = reshape(epi_params.ρˢᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "E"
        df.cases = reshape(epi_params.ρᴱᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "A"
        df.cases = reshape(epi_params.ρᴬᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "I"
        df.cases = reshape(epi_params.ρᴵᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "PH"
        df.cases = reshape(epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "PD"
        df.cases = reshape(epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "HR"
        df.cases = reshape(epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "HD"
        df.cases = reshape(epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "D"
        df.cases = reshape(epi_params.ρᴰᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    elseif compartment == "R"
        df.cases = reshape(epi_params.ρᴿᵍᵥ .* population.nᵢᵍ, G * M * T * V)
    end

    CSV.write(@sprintf("%s/output_%s_%s.csv", folder, compartment, suffix), df)
end


"""
    print_status(epi_params::Epidemic_Params,
                 population::Population_Params,
                 t::Int64)

Print the status of the epidemic spreading.
"""
function print_status(epi_params::Epidemic_Params,
                      population::Population_Params,
                      t::Int64)

    players  = sum((epi_params.ρˢᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, :] .+
                    epi_params.CHᵢᵍᵥ[:, :, t, :] ) .* population.nᵢᵍ[:, :])

    sus3 = sum((epi_params.ρˢᵍᵥ[:, :, t, 3] ) .* population.nᵢᵍ[:, :])

    infected = sum(epi_params.ρᴵᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :] .+
                   epi_params.ρᴬᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    cases3    = sum((epi_params.ρᴾᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 3]) .* population.nᵢᵍ[:, :])

    icus     = sum((epi_params.ρᴴᴿᵍᵥ[:, :, t, :] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, :]) .* population.nᵢᵍ[:, :])

    deaths   = sum(epi_params.ρᴰᵍᵥ[:, :, t, :] .* population.nᵢᵍ[:, :])

    vaccine1 = sum((epi_params.ρˢᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 1] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 1] .+
                    epi_params.CHᵢᵍᵥ[:, :, t, 1] ) .* population.nᵢᵍ[:, :]) / population.N

    vaccine2 = sum((epi_params.ρˢᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 2] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 2] .+
                    epi_params.CHᵢᵍᵥ[:, :, t, 2] ) .* population.nᵢᵍ[:, :]) / population.N

    vaccine3 = sum((epi_params.ρˢᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴱᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴬᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴵᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴾᴴᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴴᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴿᵍᵥ[:, :, t, 3] .+
                    epi_params.ρᴰᵍᵥ[:, :, t, 3] .+
                    epi_params.CHᵢᵍᵥ[:, :, t, 3] ) .* population.nᵢᵍ[:, :]) / population.N

    @printf("Time: %d, players: %.2f, sus3: %.2f, cases3: %.2f, deaths: %.2f, vaccine1 = %.2f, vaccine2: %.2f, vaccine3: %.2f\n",
            t, players, sus3, cases3, deaths, vaccine1, vaccine2, vaccine3 )
    
end