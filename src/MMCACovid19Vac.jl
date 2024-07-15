module MMCACovid19Vac

using Statistics, Dates, Printf, Logging
using CSV, NPZ, JSON, HDF5, DataStructures, DelimitedFiles, DataFrames, NetCDF


include("utils.jl")
export Epidemic_Params,
	Population_Params,
	update_epidemic_params!,
	update_population_params!,
	reset_epidemic_params!,
	reset_params!,
	set_initial_conditions!,
	correct_self_loops,
	compute_R_eff

include("mmca.jl")
export run_epidemic_spreading_mmca!,
    run_epidemic_spreading_mmca!
    

end
