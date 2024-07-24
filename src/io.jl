function validate_config(config)
    @assert !haskey(config, "epidemic_params")
    @assert !haskey(config, "population_params")
    @assert !haskey(config, "vaccination")
    @assert !haskey(config, "NPI") 
end