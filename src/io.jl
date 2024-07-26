function parse_commandline()
    s = ArgParseSettings()
    engines = ["MMCACovid19Vac"]

    @add_arg_table! s begin
        "--config", "-c"
            help = "config file (json file)"
            required = true
        "--engine", "-e"
            help = "Simulator Engine"
            range_tester = (x->x ∈ engines)
            default = "MMCACovid19Vac" 
        "--data-folder", "-d"
            help = "data folder"
            required = true
        "--instance-folder", "-i"
            help = "instance folder (experiment folder)"
            default = "." 
        "--export-compartments-full"
            help = "export compartments of simulations"
            action = :store_true
        "--export-compartments-time-t"
            help = "export compartments of simulations at a given time"
            default = nothing
            arg_type = Int
        "--initial-condition"
            help = "compartments to initialize simulation. If missing, use the seeds to initialize the simulations"
            default = ""
        "--start-date"
            help = "starting date of simulation. Overwrites the one provided in config.json"
            default = nothing
        "--end-date"
            help = "end date of simulation. Overwrites the one provided in config.json"
            default = nothing
    end
    return parse_args(s)
end

function validate_config(config, engine)
    if engine == "MMCACovid19Vac"
        @assert haskey(config, "simulation")
        @assert haskey(config, "data")
        @assert haskey(config, "epidemic_params")
        @assert haskey(config, "population_params")
        @assert haskey(config, "vaccination")
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