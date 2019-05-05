push!(LOAD_PATH, "../src/")

using NSHeat
using Logging

# Chnage here for your direcotry
ROOT_DIR = "./cool/"

function run_cool(model, core, env, var)
    @show model.modelname
        
    try
        sol = cooling(model, core, env, var)
        @show sol[end]
        if sol[end] == :Success
            write_ini(sol, model)
            output_T(sol, model, core, env, var)
            output_LC(sol, model, core, env, var)
        else
            @error "ODE solver failed" mode.modelname
        end
    catch err
        @show "Failed"
        @error "ODE solver failed" model.modelname
        @error "Error" err
    end
end
        
function main()
    println(PROGRAM_FILE," start!!")

    eos = "APR_EOS_Cat_core.dat"
    tov = "Prof_APR_Cat_1.4.dat"
    dMoverM = 0.578e-18/1.732
    del_slice = 100.0

    Tinf0 = 1.e+11
    tyr0 = 0.1
    eta_e_inf0 = 1e-30
    eta_mu_inf0 = 1e-30

    gapmodel_ns = ["a", "b", "c"]

    gapmodel_ps = ["AO", "CCDK"]

    # not used for cooling
    noneq = false
    P0 = 1e-3
    Pnow = 5.8e-3
    Pdotnow = 5.7e-20
    Znpe = 1.e-60
    Znpmu = 1.2e-60
    Znp = 4.e-61
    Wnpe = -1.5e-13
    Wnpmu = -2.e-13

    solver = "CVODE_BDF"
    tyrf = 1e7
    reltol = 1e-4
    abstol = 1e-4
    dt = 0.05

    io = open(ROOT_DIR * "log.txt", "w+")
    logger = ConsoleLogger(io)

    with_logger(logger) do
        @info "Results are placed in " ROOT_DIR

        SFtype_n = "Normal"
            
        for SFtype_p = ["Normal", "1S0"]
            if SFtype_p == "Normal"
                modelname = "$(SFtype_n)_$(SFtype_p)"
                output_dir = ROOT_DIR *modelname
                model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                              Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                              SFtype_n, "a", SFtype_p, "AO",
                                              noneq, P0, Pnow, Pdotnow,
                                              Znpe, Znpmu, Znp, Wnpe, Wnpmu,
                                              solver, tyrf, reltol, abstol, dt,
                                              output_dir)
                run_cool(model, core, env, var)
            else
                for gapmodel_p in gapmodel_ps
                    modelname = "$(SFtype_n)_$(SFtype_p)_$(gapmodel_p)"
                    output_dir = ROOT_DIR * modelname
                    model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                                  Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                                  SFtype_n, "a", SFtype_p, gapmodel_p,
                                                  noneq, P0, Pnow, Pdotnow,
                                                  Znpe, Znpmu, Znp, Wnpe, Wnpmu,
                                                  solver, tyrf, reltol, abstol, dt,
                                                  output_dir)
                    run_cool(model, core, env, var)
                end
            end

        end

        SFtype_n = "3P2m0"
        for SFtype_p = ["Normal", "1S0"]
            if SFtype_p == "Normal"
                for gapmodel_n=gapmodel_ns
                    modelname = "$(SFtype_n)_$(gapmodel_n)_$(SFtype_p)"
                    output_dir = ROOT_DIR * modelname
                    model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                                  Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                                  SFtype_n, gapmodel_n, SFtype_p, "AO",
                                                  noneq, P0, Pnow, Pdotnow,
                                                  Znpe, Znpmu, Znp, Wnpe, Wnpmu,
                                                  solver, tyrf, reltol, abstol, dt,
                                                  output_dir)
                    run_cool(model, core, env, var)
                end
            else
                for gapmodel_n=gapmodel_ns, gapmodel_p=gapmodel_ps
                    modelname = "$(SFtype_n)_$(gapmodel_n)_$(SFtype_p)_$(gapmodel_p)"
                    output_dir = ROOT_DIR * modelname
                    model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                                  Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                                  SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                                  noneq, P0, Pnow, Pdotnow,
                                                  Znpe, Znpmu, Znp, Wnpe, Wnpmu,
                                                  solver, tyrf, reltol, abstol, dt,
                                                  output_dir)
                    run_cool(model, core, env, var)
                end
            end
        end

        flush(io)
    end
    close(io)
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

