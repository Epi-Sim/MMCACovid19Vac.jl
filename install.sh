#!/bin/bash

source hpc_utils.sh

MMCACovid19Vac="https://github.com/Epi-Sim/MMCACovid19Vac.jl"

if ! in_hpc_bsc || in_hpc_wifi; then
    echo "Installing MMCACovid19Vac package..."
    if julia --project . -e "using Pkg; Pkg.add(url=\"$MMCACovid19Vac\"); Pkg.instantiate(); Pkg.precompile()"; then
        echo "Engine installed successfully."
    else
        echo "Engine installation failed. Please check the error messages above."
        exit 1
    fi
else
    echo "Skipping engine installation in this environment."
    echo "Make sure the engines are already installed !!!"
fi

echo "Compiling the package..."
if in_hpc_bsc; then
    echo "Starting interactive session"
    srun -t 00:30:00 -A $USER --qos gp_bscls -c 4 -n 1 julia install.jl
else
    julia install.jl
fi