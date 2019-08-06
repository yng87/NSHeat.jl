# NSHeat
Nuetron star cooling and heating code.

The code solves temperature evolution of spherically symmetric and isothermal NS, along with the evolution of the imbalance among chemical potentials.

The imbalance causes the internal heating because of the entropy production.

The code can also incorporate the heating caused by dark matter accretion.

# Julia installation

See the [official site](https://julialang.org/downloads/).

For mac homebrew user
```console
$ brew cask install julia
```

# Package installation

Please install ```gfortran``` in advance.

## From GitHub

```console
$ julia
julia> ]
(v1.1) pkg> add https://github.com/yng87/NSHeat.jl.git
```

Then in your julia code, you can use NSHeat by
```julia
using NSHeat
```

The dependent modules are automatically installed.

Currently, however, the module ```ODEInterface``` has a trouble finding correct path.
So please install it manually:
```console
pkg> add ODEInterface
```

When you encounter error in building NSHeat, try ```pkg> resolve```.

## Manual installation
If you build it by yourself, just download zip.
In your julia code, specify the path explicitly:
```julia
push!(LOAD_PATH, "path/to/NSHeat")
using NSHeat
```
You may need to install other packages by `add [package name]`.

# How to run

## Set Model parameters

NSHeat uses four structures to manage the calculation.

- `ModelParams`: model specification such as EOS or superfluid gap models. 

- `StarCoreParams`: T=0 quantities given by EOS such as Fermi momenta or critical temperatures of nucleons.

- `EnvelopeParams`: envelope paramters govering the surface photon emission

- `StarVariables`: temperature and chemical potentials which change by time

One can set the parameters for each struct by the following way:

1. From input cards: 
`.ini` file is supported. See `test/sample.ini` and `test/test.jl`

2. Directly in the code:
See e.g., `test/test_cool.jl`

## Solve ODE

The basic usage is read from the test codes in `test/`
(In particular, see notebook `test/test.ipynb` at first).



The calculation goes as follows:

When passing the parameters by ini file, NS parameters are set by  
```julia
model, core, env, var = setup("path/to/ini/file/")
```  
These four variables have type of sturcts defined above respectively.

Then you can solve ODE by
```julia
sol = heating(model, core, env, var)
```
or 
```julia
sol = cooling(model, core, env, var)
```

`sol` stores the solution of ODE. 

The results are placed in the directory you specified in ini card by

```julia
write_ini(sol, model) # save the model parameters
```

```julia
output_T(sol, model, core, env, var) # save temperature and chemical potential
```

```julia
output_LC(sol, model, core, env, var) # save luminosities and heat capacities
```

# ODE Solver
Julia DifferentialEquation.jl offers a lot of solvers for an ODE problem.
Among them, ```radau``` is the best for this NS evolution, in particular for late time heating problem.
