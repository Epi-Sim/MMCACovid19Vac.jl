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
	compute_R_eff

include("mmca.jl")
export run_epidemic_spreading_mmca!,
    run_epidemic_spreading_mmca!
    

include("io.jl")
export create_default_epi_params,
	create_default_population_params,
	create_default_vacparameters,
	create_default_npiparameters,
	init_pop_param_struct,
	init_epi_parameters_struct,
	init_NPI_parameters_struct,
	save_simulation_hdf5,
	save_simulation_netCDF

end
