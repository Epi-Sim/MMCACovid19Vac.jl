using Pkg

Pkg.activate(".")

Pkg.add(url="https://github.com/Epi-Sim/MMCACovid19Vac.jl")
Pkg.instantiate()
Pkg.precompile()

using PackageCompiler
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--compile", "-c"
            help = "Compile the simulator into a single precompiled excecutable"
            action = :store_true
        "--target", "-t"
            help = "Target folder where the single excecutable will be stored"
            default ="."
    end
    return parse_args(s)
end



args = parse_commandline()
if args["compile"]
    build_folder = "build"
    create_app(pwd(), build_folder, force=true)
    bin_path = abspath(joinpath(build_folder, "bin", "EpiSim"))
    symlink_path = joinpath(args["compile"], "episim")
    symlink(bin_path, symlink_path)
end