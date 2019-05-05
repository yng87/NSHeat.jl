__precompile__()

module NSHeat
"""
See test codes in test directory for how to use this module.
"""

export cooling, heating, write_ini, output_T, output_LC, setup, heating_log, heating_wimp

push!(LOAD_PATH, "./")
using DifferentialEquations
using Sundials
using LSODA
using ODEInterfaceDiffEq
using DelimitedFiles
using Dierckx

nsheat_path = joinpath(dirname(pathof(NSHeat)), "../")

include("./PhysicalConstants.jl")
include("./Simpson.jl")
include("./MyConfParser.jl")

# Basic star codes
include("./NeutronStar.jl")
include("./Setup.jl")
include("./Output.jl")
include("./SuperfluidGaps.jl")

# heat capacity
include("./SpecHeat.jl")
include("./HeatCapacity.jl")

# Neutrino emission
include("./Urca.jl")
include("./UrcaNoneq.jl")
include("./NeutrinoPBF.jl")
include("./NeutrinoLum.jl")

# Photon emission
include("./PhotonLum.jl")

# spin-down
include("./SpinDown.jl")

# ODE solver
include("./ODESolvers.jl")

# DM heating
include("./DMLum.jl")

end
