# NSHeat
Nuetron star cooling and heating code.

The code solves temperature evolution of spherically symmetric and isothermal NS, along with the evolution of the imbalance among chemical potentials.

The imbalance causes the internal heating because of the entropy production.

# Julia installation

# Package installation

# How to run

## Set Model parameters

NSHeat uses four structures to manage the calculation.

- `ModelParams`: model specification such as EOS or superfluid gap models. 

- `StarCoreParams`: T=0 quantities given by EOS such as Fermi momenta or critical temperatures of nucleons.

- `EnvelopeParams`: envelope paramters govering the surface photon emission

- `StarVariables`: temperature and chemical potentials which change by time

One can set the parameters for each struct by the following way:

1. From input cards: 
`.ini` file is supported. See `src/sample.ini` and `test/test.jl`

2. Directly in the code:
See e.g., `test/test_cool.jl`

## Solve ODE

The basic usage is read from the test codes in `test/`
(In particular, see notebook `test/test.ipynb` at first).



The calculation goes as follows:

1. When passing the parameters by ini file, NS parameters are set by

```model, core, env, var = setup(path/to/ini/file/)```.

These four variables have type of sturcts defined above respectively.

2. Then you can solve ODE by

```sol = heating(model, core, env, var)```  or ```sol = cooling(model, core, env, var)```


`sol` stores the solution of ODE. 

3. The results are placed in the directory you specified in ini card by

```write_ini(sol, model)``` (save the model parameters)

```output_T(sol, model, core, env, var)``` (save temperature and chemical potential)

```output_LC(sol, model, core, env, var)``` (save luminosities and heat capacities)

# About precision
The code works really well for cooling only problem.
Including heating triggered by non-equilibrium urca process, the differential equation becomes very stiff.
In the case of large neutron gap, the solving ODE sometimes fails or continues to run forever.
Then one must tune the tolerance giving to the ODE solver.
