push!(LOAD_PATH, "../src/")

using NSHeat
using Logging

# chnage hear to your directory
ROOT_DIR = "./heat/"

function run_heat(model, core, env, var)
    @show model.modelname
        
    try
        sol = heating(model, core, env, var)
        @show sol[end]
        if sol[end] == :Success
            write_ini(sol, model)
            output_T(sol, model, core, env, var)
            output_LC(sol, model, core, env, var)
        else
            write_ini(sol, model)
            @error "ODE solver failed" model.modelname
            @error "ODE" sol[end]
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

    # Millisecond pulsarp
    eos = "../EOS_data/APR_EOS_Cat_core.dat"
    del_slice = 10.0

    Tinf0 = 1.e+10
    tyr0 = 0.1
    eta_e_inf0 = 1e-30
    eta_mu_inf0 = 1e-30
    SFtype_n = "3P2m0"
    gapmodel_ns = ["a", "b", "c"]
    #gapmodel_ns = ["c"]
    SFtype_p = "1S0"
    gapmodel_ps = ["AO", "CCDK"]
    
    # not used for cooling
    noneq = true
    Znpe = Dict("1.4"=>1e-60, "1.8"=>6e-61)
    Znpmu = Dict("1.4"=>1.2e-60, "1.8"=>7e-61)
    Znp = Dict("1.4"=>4e-61, "1.8"=>2e-61)
    Wnpe = Dict("1.4"=>-1.5e-13, "1.8"=>-1.4e-13)
    Wnpmu = Dict("1.4"=>-2e-13, "1.8"=>-1.8e-13)

    alpha = 10.0
    beta = 100.0
    dt = 0.05
    solver = "CVODE_BDF"

    io = open(ROOT_DIR * "log.txt", "w+")
    logger = ConsoleLogger(io)

    with_logger(logger) do
        @info "Results are placed in " ROOT_DIR

        #For classical pulsar
        masses = ["1.4", "1.8"]
        dMs = [1e-7, 1e-15]
        Pnow = 1.0
        Pdotnow = 1e-15
        P0s = [1e-3, 1e-2, 1e-1]
        tyrf = 1e9
        reltol=1e-3
        abstol=1e-6
        for gapmodel_n=gapmodel_ns, gapmodel_p=gapmodel_ps, mass=masses, dMoverM=dMs, P0=P0s
            tov = "../TOV_data/Profile/Prof_APR_Cat_$(mass).dat"
            modelname = "CP_$(gapmodel_n)_$(gapmodel_p)_$(mass)_$(dMoverM)_$(P0)"
            output_dir = ROOT_DIR * modelname
            model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                          Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                          SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                          noneq, P0, Pnow, Pdotnow,
                                          Znpe[mass], Znpmu[mass], Znp[mass], Wnpe[mass], Wnpmu[mass],
                                          solver, tyrf, reltol, abstol, dt,
                                          alpha, beta,
                                          output_dir)
            run_heat(model, core, env, var)
        end

        # # MSP
        mass = "1.4"
        Pnow = 5.8e-3
        Pdotnow = 5.7e-20
        P0s = [1e-3, 0.5e-3]
        tyrf = 1e10
        reltol = 1e-4
        abstol = 1e-4
        dMoverM = 1e-7
        for gapmodel_n=gapmodel_ns, gapmodel_p=gapmodel_ps, P0=P0s
            tov = "../TOV_data/Profile/Prof_APR_Cat_1.4.dat"
            modelname = "MSP_$(gapmodel_n)_$(gapmodel_p)_$(P0)"
            output_dir = ROOT_DIR * modelname
            model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                          Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                          SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                          noneq, P0, Pnow, Pdotnow,
                                          Znpe[mass], Znpmu[mass], Znp[mass], Wnpe[mass], Wnpmu[mass],
                                          solver, tyrf, reltol, abstol, dt,
                                          alpha, beta,
                                          output_dir)
            run_heat(model, core, env, var)
        end
        flush(io)
    end
    close(io)
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

