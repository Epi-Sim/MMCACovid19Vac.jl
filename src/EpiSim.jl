module EpiSim

using MMCACovid19Vac
using Dates, Logging
using HDF5, DataFrames, NetCDF


include("io.jl")
export create_default_epi_params,
    create_default_population_params,
    create_default_vacparameters,
    create_default_npiparameters,
    update_config!,
    init_pop_param_struct,
    init_epi_parameters_struct,
    init_NPI_parameters_struct,
    save_simulation_hdf5,
    save_simulation_netCDF
end # module EpiSim
