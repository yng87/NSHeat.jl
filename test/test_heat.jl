push!(LOAD_PATH, "../src/")

using Setup
using ODESolvers
using Output
using Logging

function run_heat(model, core, env, var, reltol, abstol)
    @show model.modelname
        
    try
        sol = heating(model, core, env, var, reltol, abstol)
        @show sol.retcode
        if sol.retcode == :Success
            write_ini(sol, model)
            output_T(sol, model, core, env, var)
            output_LC(sol, model, core, env, var)
        else
            write_ini(sol, model)
            @error "ODE solver failed" model.modelname
            @error "ODE" sol.retcode
        end
    catch err
        #write_ini(sol, model)
        @show "Failed"
        @error "ODE solver failed" model.modelname
        @error "Error" err
    end
end
        
function main()
    println(PROGRAM_FILE," start!!")

    # Millisecond pulsar
    eos = "../EOS_data/APR_EOS_Cat_core.dat"
    tov = "../TOV_data/Profile/Prof_APR_Cat_1.4.dat"
    dMoverM = 1e-7
    del_slice = 100.0

    Tinf0 = 1.e+11
    tyr0 = 0.1
    eta_e_inf0 = 1e-30
    eta_mu_inf0 = 1e-30
    SFtype_n = "3P2m0"
    gapmodel_ns = ["a", "b", "c"]
    SFtype_p = "1S0"
    gapmodel_ps = ["AO", "CCDK"]

    # not used for cooling
    noneq = true
    Pnow = 5.8e-3
    Pdotnow = 5.7e-20
    Znpe = Dict("1.4"=>1e-60, "1.8"=>6e-61)
    Znpmu = Dict("1.4"=>1.2e-60, "1.8"=>7e-61)
    Znp = Dict("1.4"=>4e-61, "1.8"=>2e-61)
    Wnpe = Dict("1.4"=>-1.5e-13, "1.8"=>-1.4e-13)
    Wnpmu = Dict("1.4"=>-2e-13, "1.8"=>-1.8e-13)
    @show Znpe
    solver = "CVODE_BDF"
    tyrf = 1e10

    ROOT_DIR = homedir() * "/Dropbox/MyWorks/rotochemical/NSHeat/heating_nonzeroT/"

    io = open(ROOT_DIR * "log.txt", "w+")
    logger = ConsoleLogger(io)

    with_logger(logger) do
        @info "Results are placed in " ROOT_DIR

        # For classical pulsar
        masses = ["1.4", "1.8"]
        dMs = [1e-7, 1e-15]
        Pnow = 1.0
        Pdotnow = 1e-15
        P0s = [1e-3, 1e-2, 1e-1]

        for gapmodel_n=gapmodel_ns, gapmodel_p=gapmodel_ps, mass=masses, dMoverM=dMs, P0=P0s
            tov = "../TOV_data/Profile/Prof_APR_Cat_$(mass).dat"
            modelname = "CP_$(gapmodel_n)_$(gapmodel_p)_$(mass)_$(dMoverM)_$(P0)"
            output_dir = ROOT_DIR * modelname
            model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                          Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                          SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                          noneq, P0, Pnow, Pdotnow,
                                          Znpe[mass], Znpmu[mass], Znp[mass], Wnpe[mass], Wnpmu[mass],
                                          solver, tyrf,
                                          output_dir)
            run_heat(model, core, env, var, 1e-4, 1e-4)
        end

        # MSP
        mass = "1.4"
        P0s = [1e-3, 0.5e-3]
        for gapmodel_n=gapmodel_ns, gapmodel_p=gapmodel_ps, P0=P0s
            modelname = "MSP_$(gapmodel_n)_$(gapmodel_p)_$(P0)"
            output_dir = ROOT_DIR * modelname
            model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                          Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                          SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                          noneq, P0, Pnow, Pdotnow,
                                          Znpe[mass], Znpmu[mass], Znp[mass], Wnpe[mass], Wnpmu[mass],
                                          solver, tyrf,
                                          output_dir)
            run_heat(model, core, env, var, 3e-3, 3e-3)
        end
        flush(io)
    end
    close(io)
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

