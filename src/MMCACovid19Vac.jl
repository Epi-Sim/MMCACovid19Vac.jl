module MMCACovid19Vac

using Statistics, Dates, Printf, Logging, LinearAlgebra
using CSV, NPZ, JSON, HDF5, DataStructures, DelimitedFiles, DataFrames, NetCDF


include("utils.jl")
export Epidemic_Params,
	Population_Params,
	NPI_Params,
	init_pop_param_struct,
	init_epi_parameters_struct,
	init_NPI_parameters_struct,
	update_epidemic_params!,
	reset_epidemic_compartments!,
	reset_params!,
	set_initial_conditions!,
	set_compartments!,
	correct_self_loops


include("mmca.jl")
export run_epidemic_spreading_mmca!,
	update_population_params!,
	run_epidemic_spreading_mmca!


include("io.jl")
export create_default_epi_params,
	create_default_population_params,
	create_default_vacparameters,
	create_default_npiparameters,
	create_config_template,
	save_simulation_hdf5,
	save_simulation_netCDF,
	store_R_eff,
	print_status


include("compute.jl")
export compute_R_eff,
	optimal_vaccination_distribution


end
