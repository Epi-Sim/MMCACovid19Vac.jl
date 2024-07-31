function validate_config(config, engine)
    if engine == "MMCACovid19Vac"
        @assert haskey(config, "simulation")
        @assert haskey(config, "data")
        @assert haskey(config, "epidemic_params")
        @assert haskey(config, "population_params")
        @assert haskey(config, "vaccination")
        @assert haskey(config, "NPI")
    elseif engine == "MMCACovid19"
        @assert haskey(config, "simulation")
        @assert haskey(config, "data")
        @assert haskey(config, "epidemic_params")
        @assert haskey(config, "population_params")
        @assert haskey(config, "NPI") 
    end
end


function update_config!(config, cmd_line_args)
    # Define dictionary containing epidemic parameters

    # overwrite config with command line
    if cmd_line_args["start-date"] !== nothing
        config["simulation"]["start_date"] = cmd_line_args["start-date"]
    end
    if cmd_line_args["end-date"] !== nothing
        config["simulation"]["end_date"] = cmd_line_args["end-date"]
    end
    if cmd_line_args["export-compartments-time-t"] !== nothing
        config["simulation"]["export_compartments_time_t"] = cmd_line_args["export-compartments-time-t"]
    end
    if cmd_line_args["export-compartments-full"] == true
        config["simulation"]["export_compartments_full"] = true
    end

    nothing
end