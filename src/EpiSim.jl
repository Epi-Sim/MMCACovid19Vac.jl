module EpiSim

using MMCACovid19Vac
using ArgParse
using Dates, Logging, Printf
using HDF5, DataFrames, NetCDF

import JSON
import CSV


include("commands.jl")
include("engine.jl")
include("io.jl")

const ENGINES  = ["MMCACovid19Vac"]
const COMMANDS = ["run", "setup", "init"]

function julia_main()::Cint
    try
        args = parse_command_line()
        command = args["%COMMAND%"]
        engine = args["engine"]
    
        # Check if the provided command is in the list of accepted commands
        if !(command in COMMANDS)
            println("Unknown command: $command")
            println("Accepted commands are: $(join(COMMANDS, ", "))")
            return 1
        end

        if command == "run"
            execute_run(args["run"], engine)
        elseif command == "setup"
            execute_setup(args["setup"], engine)
        elseif command == "init"
            execute_init(args)
        end

    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

end # module EpiSim
