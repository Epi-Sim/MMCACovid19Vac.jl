abstract type AbstractEngine end
abstract type AbstractOutputFormat end

# Add to this as we add more engines
struct MMCACovid19VacEngine <: AbstractEngine end
struct MMCACovid19Engine <: AbstractEngine end

struct NetCDFFormat <: AbstractOutputFormat end
struct HDF5Format <: AbstractOutputFormat end