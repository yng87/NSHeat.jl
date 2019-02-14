__precompile__()

module NSHeat

export cooling, heating, write_ini, output_T, output_LC, setup

push!(LOAD_PATH, "./")
using DifferentialEquations
using Sundials
using LSODA
using ODEInterfaceDiffEq
using DelimitedFiles
using Dierckx
using ConfParser

include("./PhysicalConstants.jl")
include("./Simpson.jl")

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

end
