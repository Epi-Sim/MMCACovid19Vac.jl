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

    function create_epi_params_dict(G::Int)
        epiparams_dict = Dict()
        epiparams_dict["scale_β"] = 0.5
        epiparams_dict["βᴬ"] = 0.05
        epiparams_dict["βᴵ"] = 0.09
        epiparams_dict["ηᵍ"] = ones(G) * 0.275
        epiparams_dict["αᵍ"] = ones(G) * 0.65
        epiparams_dict["μᵍ"] = ones(G) * 0.3
        epiparams_dict["θᵍ"] = zeros(G)
        epiparams_dict["γᵍ"] = ones(G) * 0.03
        epiparams_dict["ζᵍ"] = ones(G) * 0.12
        epiparams_dict["λᵍ"] = ones(G) * 0.275
        epiparams_dict["ωᵍ"] = ones(G) * 0.1
        epiparams_dict["ψᵍ"] = ones(G) * 0.14
        epiparams_dict["χᵍ"] = ones(G) * 0.047
        epiparams_dict["Λ"] = 0.02
        epiparams_dict["Γ"] = 0.01
        epiparams_dict["rᵥ"] = ones(G) * 0.6
        epiparams_dict["kᵥ"] = ones(G) * 0.4
        epiparams_dict["risk_reduction_dd"] = 0.0
        epiparams_dict["risk_reduction_h"] = 0.1
        epiparams_dict["risk_reduction_d"] = 0.05
        return epiparams_dict
    end

    function create_default_population_params()
        population = Dict()
        population["G_labels"] = ["Y", "M", "O"]
        population["C"] = [ 0.598  0.38486 0.01714 ;
                            0.244  0.721   0.0353;
                            0.1919 0.5705  0.2376
                        ]
        population["kᵍ"] = [11.8, 13.3, 6.76]
        population["kᵍ_h"] = [3.15, 3.17, 3.28]
        population["kᵍ_w"] = [1.72, 5.18, 0.0]
        population["pᵍ"] = [0.0, 1.0, 0.00]
        population["ξ"] = 0.01
        population["σ"] = 2.5

        # Mobility variation parameters
        population["mobility_variation_type"] = "none"
        population["mobility_sigma_O"] = 0.0
        population["mobility_sigma_D"] = 0.0
        
        return population
    end

    function create_population_params_dict(G::Int)
        population = Dict()
        population["G_labels"] = ["G$(i)" for i in range(1,G)]
        population["C"] = C = Matrix{Float64}(I, G, G)
        population["kᵍ"] = ones(G) * 12
        population["kᵍ_h"] = ones(G) * 3.15
        population["kᵍ_w"] = ones(G) * 2
        population["pᵍ"] = ones(G)
        population["ξ"] = 0.01
        population["σ"] = 2.5

        # Mobility variation parameters
        population["mobility_variation_type"] = "none"  # "none", "ipfp_simple", "external_netcdf"
        population["mobility_sigma_O"] = 0.0
        population["mobility_sigma_D"] = 0.0

        # External mobility NetCDF parameters
        population["mobility_external_file"] = nothing
        population["mobility_external_variable"] = "mobility"
        population["mobility_validate"] = true
        population["mobility_fill_missing"] = false
        population["mobility_time_alignment"] = "start"  # "start", "end", "nearest"

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

    function create_vacparameters_dict(G::Int)
        vacparams_dict = Dict()
        vacparams_dict["ϵᵍ"] = ones(G) * 0.3
        vacparams_dict["percentage_of_vacc_per_day"] = 0.005
        vacparams_dict["start_vacc"] = 1
        vacparams_dict["dur_vacc"] = 5
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

    function create_config_template(G::Int)
        # Define dictionary containing epidemic parameters
        config = Dict()
        config["epidemic_params"] = create_epi_params_dict(G)
        config["population_params"] = create_population_params_dict(G)
        config["vaccination"] = create_vacparameters_dict(G)
        config["NPI"] = create_default_npiparameters()
        return config
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
        compartments[:, :, :, :, 11] .= epi_params.CHᵢᵍᵥ .* population.nᵢᵍ
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
                                        G_coords= String[],
                                        M_coords = String[],
                                        T_coords = String[]
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
                                G_coords=String[], M_coords=String[], T_coords=String[]
                                )
    G = population.G
    M = population.M
    T = epi_params.T
    V = epi_params.V

    if length(G_coords) != G
        G_coords = collect(1:G)
    end
    if length(M_coords) != M
        M_coords = collect(1:M)
    end
    if length(T_coords) != T
        T_coords = collect(1:T) 
    end
    
    V_coords = epi_params.VaccLabels

    # Prepare data dictionary
    data_dict = Dict()
    data_dict["S"]  = epi_params.ρˢᵍᵥ .* population.nᵢᵍ
    data_dict["E"]  = epi_params.ρᴱᵍᵥ  .* population.nᵢᵍ
    data_dict["A"]  = epi_params.ρᴬᵍᵥ  .* population.nᵢᵍ
    data_dict["I"]  = epi_params.ρᴵᵍᵥ  .* population.nᵢᵍ
    data_dict["PH"] = epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ
    data_dict["PD"] = epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ
    data_dict["HR"] = epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ
    data_dict["HD"] = epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ
    data_dict["R"]  = epi_params.ρᴿᵍᵥ  .* population.nᵢᵍ
    data_dict["D"]  = epi_params.ρᴰᵍᵥ  .* population.nᵢᵍ
    data_dict["CH"] = epi_params.CHᵢᵍᵥ .* population.nᵢᵍ

    # Remove file if exists
    isfile(output_fname) && rm(output_fname)

    # Create dataset using NCDatasets
    ds = NCDataset(output_fname, "c")

    try
        # Define dimensions
        defDim(ds, "G", G)
        defDim(ds, "M", M)
        defDim(ds, "T", T)
        defDim(ds, "V", V)

        # Define variables
        defVar(ds, "S", Float64, ("G", "M", "T", "V"))
        defVar(ds, "E", Float64, ("G", "M", "T", "V"))
        defVar(ds, "A", Float64, ("G", "M", "T", "V"))
        defVar(ds, "I", Float64, ("G", "M", "T", "V"))
        defVar(ds, "PH", Float64, ("G", "M", "T", "V"))
        defVar(ds, "PD", Float64, ("G", "M", "T", "V"))
        defVar(ds, "HR", Float64, ("G", "M", "T", "V"))
        defVar(ds, "HD", Float64, ("G", "M", "T", "V"))
        defVar(ds, "R", Float64, ("G", "M", "T", "V"))
        defVar(ds, "D", Float64, ("G", "M", "T", "V"))
        defVar(ds, "CH", Float64, ("G", "M", "T", "V"))

        # Set attributes
        ds["S"].attrib["description"] = "Suceptibles"
        ds["E"].attrib["description"] = "Exposed"
        ds["A"].attrib["description"] = "Asymptomatic"
        ds["I"].attrib["description"] = "Infected"
        ds["PH"].attrib["description"] = "Pre-hospitalized"
        ds["PD"].attrib["description"] = "Pre-deceased"
        ds["HR"].attrib["description"] = "Hospitalized-good"
        ds["HD"].attrib["description"] = "Hospitalized-bad"
        ds["R"].attrib["description"] = "Recovered"
        ds["D"].attrib["description"] = "Dead"
        ds["CH"].attrib["description"] = "Confined"

        # Write data
        for var_label in epi_params.CompLabels
            ds[var_label][:,:,:,:] = data_dict[var_label]
        end
    finally
        close(ds)
    end

end

"""
save_observables_netCDF(epi_params::Epidemic_Params, 
                                    population::Population_Params,
                                    output_fname::String;
                                    G_coords= String[],
                                    M_coords = String[],
                                    T_coords = String[]
                                    )

    Calculate and store observable varialbes from a simulation.
    Observables include new infected, new hospitalizations and new death.
    Vaccination state is ignored by summing over vaccination state V.
    - New daily infections         :   ρᴬᵍ .* nᵢᵍ .* αᵍ
    - New daily hospitalizationsare:   ρᴬᵍ .* nᵢᵍ .* μᵍ .* (1 .- θᵍ) .* γᵍ
    - New daily death              :   diff(ρᴰᵍ)

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
function save_observables_netCDF(   epi_params::Epidemic_Params, 
                                    population::Population_Params,
                                    output_fname::String;
                                    G_coords=String[], M_coords=String[], T_coords=String[]
                                )
    G = population.G
    M = population.M
    T = epi_params.T

    if length(G_coords) != G
        G_coords = collect(1:G)
    end
    if length(M_coords) != M
        M_coords = collect(1:M)
    end
    if length(T_coords) != T
        T_coords = collect(1:T) 
    end

    # Prepare data dictionary
    data_dict = Dict()
    data_dict["new_infected"] = sum((epi_params.ρᴬᵍᵥ  .* population.nᵢᵍ), dims=(4))[:,:,:,1] .* epi_params.αᵍ

    hosp_rates = epi_params.μᵍ .* (1 .- epi_params.θᵍ) .* epi_params.γᵍ
    hosp_rates = reshape(hosp_rates, 3, 1, 1, 3)

    data_dict["new_hospitalized"] = sum(  ((epi_params.ρᴵᵍᵥ  .* population.nᵢᵍ) .* hosp_rates), dims=4)[:,:,:,1]

    D = sum(epi_params.ρᴰᵍᵥ, dims=4)[:,:,:,1]
    data_dict["new_deaths"] = zeros(size(D))
    data_dict["new_deaths"][:, :, 2:end] = diff((D .* population.nᵢᵍ), dims=3)

    τ = 14

    # Remove file if exists
    isfile(output_fname) && rm(output_fname)

    # Create dataset using NCDatasets
    ds = NCDataset(output_fname, "c")

    try
        # Define dimensions
        defDim(ds, "G", G)
        defDim(ds, "M", M)
        defDim(ds, "T", T)

        # Define variables
        defVar(ds, "new_infected", Float64, ("G", "M", "T"))
        defVar(ds, "new_hospitalized", Float64, ("G", "M", "T"))
        defVar(ds, "new_deaths", Float64, ("G", "M", "T"))

        # Set attributes
        ds["new_infected"].attrib["description"] = "Daily infections"
        ds["new_hospitalized"].attrib["description"] = "Daily hospitalizations"
        ds["new_deaths"].attrib["description"] = "Daily deaths"

        # Write data
        ds["new_infected"][:,:,:] = data_dict["new_infected"]
        ds["new_hospitalized"][:,:,:] = data_dict["new_hospitalized"]
        ds["new_deaths"][:,:,:] = data_dict["new_deaths"]

        # Handle R_eff if T > τ
        if T <= τ
            @warn "The number of time steps T ($(T)) is smaller than τ ($(τ)), R is not computable"
        else
            data_dict["R_eff"] = zeros(Float64, G, M, T)
            Rᵢᵍ_eff = compute_R_eff_matrix(epi_params, population, τ)[6]
            data_dict["R_eff"][:, :, τ+1:end] .= Rᵢᵍ_eff

            defVar(ds, "R_eff", Float64, ("G", "M", "T"))
            ds["R_eff"].attrib["description"] = "Effective reproduction number over an infectious period of 14 days"
            ds["R_eff"][:,:,:] = data_dict["R_eff"]
        end
    finally
        close(ds)
    end
end


### ----------------------------------------------------------------------------
### MOBILITY NETCDF I/O FUNCTIONS
### ----------------------------------------------------------------------------

"""
    parse_date_strings(date_strs::Array{String, 1})

Parse date strings from NetCDF coordinate to Date objects.

# Arguments
- `date_strs::Array{String, 1}`: Array of date strings (e.g., "2020-03-10")

# Returns
- `dates::Array{Date, 1}`: Array of Date objects
"""
function parse_date_strings(date_strs::Array{String, 1})
    return [Date(d) for d in date_strs]
end


"""
    align_dates(netcdf_dates::Array{Date, 1},
               sim_start_date::Date,
               sim_end_date::Date;
               alignment::String="start",
               fill_missing::Bool=false)

Find indices in NetCDF dates that match the simulation period.

# Arguments
- `netcdf_dates::Array{Date, 1}`: Dates from NetCDF file
- `sim_start_date::Date`: Simulation start date
- `sim_end_date::Date`: Simulation end date
- `alignment::String="start"`: Alignment mode - "start", "end", or "nearest"
- `fill_missing::Bool=false`: Whether to interpolate missing dates

# Returns
- `indices::Array{Int, 1}`: Indices of NetCDF dates to use
- `aligned_dates::Array{Date, 1}`: Corresponding dates
"""
function align_dates(netcdf_dates::Array{Date, 1},
                    sim_start_date::Date,
                    sim_end_date::Date;
                    alignment::String="start",
                    fill_missing::Bool=false)
    T = (sim_end_date - sim_start_date).value + 1
    indices = zeros(Int, T)
    aligned_dates = Array{Date, 1}(undef, T)

    for (t_idx, sim_date) in enumerate(sim_start_date:Day(1):sim_end_date)
        # Find matching or nearest date in NetCDF
        if fill_missing
            # Find nearest date
            diffs = [abs((d - sim_date).value) for d in netcdf_dates]
            nearest_idx = argmin(diffs)
            indices[t_idx] = nearest_idx
            aligned_dates[t_idx] = netcdf_dates[nearest_idx]
        else
            # Find exact match
            exact_idx = findfirst(==(sim_date), netcdf_dates)
            if exact_idx === nothing
                error("Date $(sim_date) not found in NetCDF file. Use fill_missing=true to interpolate.")
            end
            indices[t_idx] = exact_idx
            aligned_dates[t_idx] = netcdf_dates[exact_idx]
        end
    end

    return indices, aligned_dates
end


"""
    load_mobility_netcdf(filename::String,
                         population::Population_Params,
                         sim_start_date::Date,
                         sim_end_date::Date;
                         varname::String="mobility",
                         validate::Bool=true,
                         fill_missing::Bool=false,
                         alignment::String="start",
                         verbose::Bool=true)

Load mobility time series from NetCDF file and convert to sparse edgelist format.

# Expected NetCDF Format
```
Dimensions:
  date: T (timesteps)
  origin: M (patches)
  destination: M (patches)

Variables:
  mobility(date, origin, destination): Float32/Float64

Coordinates:
  date: date strings (e.g., "2020-03-10")
  origin: 1 to M
  destination: 1 to M
```

# Arguments
- `filename::String`: Path to NetCDF file
- `population::Population_Params`: Population params struct with edgelist
- `sim_start_date::Date`: Simulation start date
- `sim_end_date::Date`: Simulation end date
- `varname::String="mobility"`: Name of mobility variable in NetCDF
- `validate::Bool=true`: Whether to validate the loaded series
- `fill_missing::Bool=false`: Whether to interpolate missing dates
- `alignment::String="start"`: Date alignment mode
- `verbose::Bool=true`: Whether to print progress

# Returns
- `R_series::Array{Float64, 2}`: Mobility series of shape (T, L) where L is edges
- `dates::Array{Date, 1}`: Aligned dates for each timestep
"""
function load_mobility_netcdf(filename::String,
                             population::Population_Params,
                             sim_start_date::Date,
                             sim_end_date::Date;
                             varname::String="mobility",
                             validate::Bool=true,
                             fill_missing::Bool=false,
                             alignment::String="start",
                             verbose::Bool=true)

    if verbose
        @info "Loading mobility from NetCDF: $(filename)"
    end

    M = population.M
    L = length(population.edgelist[:, 1])

    # Open NetCDF file using NCDatasets (wraps HDF5, better ARM64 support)
    ds = NCDataset(filename)

    try
        # Read mobility variable
        if !haskey(ds, varname)
            error("Variable '$(varname)' not found in NetCDF file")
        end

        mobility_dense = ds[varname][:,:,:]
        date_strings = ds["date"][:]

        # Get dimensions from the data array
        T_netcdf, M_netcdf, M_dest = size(mobility_dense)

        if verbose
            @info "  NetCDF dimensions: date=$(T_netcdf), origin=$(M_netcdf), destination=$(M_dest)"
        end

        # Validate M matches
        if M_netcdf != M
            error("Dimension mismatch: NetCDF has M=$(M_netcdf), but population has M=$(M)")
        end

        # Parse date strings to Date objects
        netcdf_dates = parse_date_strings(String.(date_strings))

        # Align dates with simulation period
        T_sim = (sim_end_date - sim_start_date).value + 1
        indices, aligned_dates = align_dates(netcdf_dates, sim_start_date, sim_end_date,
                                              alignment=alignment, fill_missing=fill_missing)

        if verbose
            @info "  Aligned $(length(indices)) timesteps from $(aligned_dates[1]) to $(aligned_dates[end])"
        end

        # Validate that external mobility covers the full simulation period
        T_loaded = length(indices)
        if T_loaded < T_sim
            loaded_start = aligned_dates[1]
            loaded_end = aligned_dates[end]
            error("External mobility data insufficient for full simulation. " *
                  "Loaded $(T_loaded) timesteps ($(loaded_start) to $(loaded_end)), " *
                  "but simulation requires $(T_sim) timesteps ($(sim_start_date) to $(sim_end_date)). " *
                  "Provide mobility data covering the entire simulation period.")
        end

        # Convert dense (T, M, M) to sparse (T, L) using edgelist
        R_series = zeros(Float64, T_sim, L)
        edgelist = population.edgelist

        for t in 1:T_sim
            nc_idx = indices[t]
            for e in 1:L
                origin = edgelist[e, 1]
                dest = edgelist[e, 2]
                R_series[t, e] = mobility_dense[nc_idx, origin, dest]
            end
        end

        # Validate if requested
        if validate
            is_valid, errors = validate_mobility_series(R_series, population.edgelist,
                                                        M, T_sim, verbose=verbose)
            if !is_valid
                @error "Mobility validation failed"
                for err in errors
                    @error "  - $(err)"
                end
                error("Invalid mobility series")
            end
        end

        if verbose
            @info "  Loaded mobility series: shape=$(size(R_series)), edges=$(L)"
        end

        return R_series, aligned_dates

    finally
        close(ds)
    end
end


"""
    save_mobility_netcdf(R_series::Array{Float64, 2},
                         edgelist::Array{Int64, 2},
                         dates::Array{Date, 1},
                         filename::String;
                         varname::String="mobility",
                         compress::Int=4)

Save mobility time series to NetCDF file.

Converts sparse (T, L) format to dense (T, M, M) and saves to NetCDF.

# Arguments
- `R_series::Array{Float64, 2}`: Mobility series of shape (T, L)
- `edgelist::Array{Int64, 2}`: Sparse edgelist of shape (L, 2)
- `dates::Array{Date, 1}`: Dates for each timestep
- `filename::String`: Output NetCDF filename
- `varname::String="mobility"`: Name of mobility variable
- `compress::Int=4`: Compression level (0-9, default 4). Set to 0 to disable.

Note: Uses NCDatasets.jl which wraps HDF5 for better ARM64 (Apple Silicon) support.
"""
function save_mobility_netcdf(R_series::Array{Float64, 2},
                             edgelist::Array{Int64, 2},
                             dates::Array{Date, 1},
                             filename::String;
                             varname::String="mobility",
                             compress::Union{Bool, Int}=4)

    T, L = size(R_series)
    M = maximum(edgelist)

    @info "Saving mobility to NetCDF: $(filename)"
    @info "  Dimensions: T=$(T), M=$(M), L=$(L)"

    # Convert sparse (T, L) to dense (T, M, M)
    mobility_dense = zeros(Float64, T, M, M)
    for t in 1:T
        for e in 1:L
            origin = edgelist[e, 1]
            dest = edgelist[e, 2]
            mobility_dense[t, origin, dest] = R_series[t, e]
        end
    end

    # Create NetCDF file using NCDatasets (wraps HDF5, better ARM64 support)
    ds = NCDataset(filename, "c"; attrib=Dict("description" => "Mobility time series for MMCA model"))

    try
        # Define dimensions
        defDim(ds, "date", T)
        defDim(ds, "origin", M)
        defDim(ds, "destination", M)

        # Define and write date variable with string format
        date_strs = [Dates.format(d, "yyyy-mm-dd") for d in dates]
        defVar(ds, "date", String, ("date",))
        ds["date"][:] = date_strs

        # Define and write mobility variable
        defVar(ds, varname, Float64, ("date", "origin", "destination"))
        ds[varname].attrib["description"] = "Mobility flow probability"
        ds[varname].attrib["units"] = "probability"
        ds[varname][:,:,:] = mobility_dense
    finally
        close(ds)
    end

    @info "Mobility saved to $(filename)"
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

        formatted_output = @sprintf("Time: %d, players: %.2f, sus3: %.2f, cases3: %.2f, deaths: %.2f, vaccine1 = %.2f, vaccine2: %.2f, vaccine3: %.2f\n",
        t, players, sus3, cases3, deaths, vaccine1, vaccine2, vaccine3 )
        @info "$(formatted_output) "    
    end
