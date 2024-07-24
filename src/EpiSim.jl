module EpiSim

using MMCACovid19Vac
using ArgParse
import JSON
import CSV
using Dates, Logging
using HDF5, DataFrames, NetCDF

include("engine.jl")

function julia_main()::Cint
    try
        args = parse_commandline()
        engine = args["engine"]
        
        if engine == "MMCACovid19Vac"
            run_MMCACovid19Vac(args)
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
