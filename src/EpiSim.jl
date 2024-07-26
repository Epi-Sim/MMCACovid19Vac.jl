module EpiSim

using MMCACovid19Vac
using ArgParse
import JSON
import CSV
using Dates, Logging
using HDF5, DataFrames, NetCDF

include("engine.jl")
include("io.jl")

function julia_main()::Cint
    try
        args = parse_commandline()
        
        engine        = args["engine"]
        data_path     = args["data-folder"]
        config_fname  = args["config"]
        instance_path = args["instance-folder"]
        init_condition_path = args["initial-condition"]
        
        config = JSON.parsefile(config_fname);
        update_config!(config, args)

        @assert isfile(config_fname);
        @assert isdir(data_path);
        @assert isdir(instance_path);
        
        validate_config(config, engine)

        if engine == "MMCACovid19Vac"
            run_MMCACovid19Vac(config, data_path, instance_path, init_condition_path)
        else
            println("Error unknown engine $(engine)")
            return 1
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

end # module EpiSim
