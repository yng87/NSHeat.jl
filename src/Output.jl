module Output

export output_T, output_LC, write_ini

push!(LOAD_PATH, "./")

using DifferentialEquations
using NeutronStar
using HeatCapacity
using NeutrinoLum
using PhotonLum
using SuperfluidGaps
using SpinDown
using ConfParser

function output_T(sol, model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    if isdir(model.output_dir) == false
        mkdir(model.output_dir)
    end
    
    filepath = model.output_dir * "/temperature.dat"

    open(filepath, "w") do io
        println(io, "# t[yr] Teff_inf[K] Tinf[K] eta_e_inf[erg] eta_mu_inf[erg]")
        for i in 1:length(sol[1])
            if model.noneq == true
                var.t = sol[1][i] #yr
                var.Tinf = sol[2][i]
                var.eta_e_inf = sol[3][i] #erg
                var.eta_mu_inf = sol[4][i] #erg
                Teff_inf = get_Teff_inf(model, env, var)
                println(io, "$(var.t) $(Teff_inf) $(var.Tinf) $(var.eta_e_inf) $(var.eta_mu_inf)")
            else
                var.Tinf = sol[2][i]
                Teff_inf = get_Teff_inf(model, env, var)
                println(io, "$(sol[1][i]) $(Teff_inf) $(var.Tinf) $(var.eta_e_inf) $(var.eta_mu_inf)")
            end

        end
    end
end

function output_LC(sol, model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    if isdir(model.output_dir) == false
        mkdir(model.output_dir)
    end
    
    filepath_L = model.output_dir * "/luminosity.dat"
    filepath_C = model.output_dir * "/capacity.dat"

    ioL = open(filepath_L, "w")
    ioC = open(filepath_C, "w")

    println(ioL, "# t[yr] Urca PBF Photon Heat")
    println(ioL, "# Luminosity unit = [erg/s]")
    println(ioC, "# t[yr] n p e mu")
    println(ioC, "# Heat Capacity unit = [erg/K]")
    for i in 1:length(sol[1])
        if model.noneq == true
            var.t = sol[1][i] #yr
            var.Tinf = sol[2][i]
            var.eta_e_inf = sol[3][i] #erg
            var.eta_mu_inf = sol[4][i] #erg
        else
            var.t = sol[1][i] #yr
            var.Tinf = sol[2][i]
        end
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)
        set_Omega(model, var)
        set_Omega_dot(model, var)
        
        Lurca = L_murca_n_e(model, core, var, model.noneq) + L_murca_n_mu(model, core, var, model.noneq) + L_murca_p_e(model, core, var, model.noneq) + L_murca_p_mu(model, core, var, model.noneq)
        
        LPBF = 0.0
        if lowercase(model.SFtype_n) != "normal"
            LPBF += L_PBF_n(model, core, var)
        end
        if lowercase(model.SFtype_p) != "normal"
            LPBF += L_PBF_p(model, core, var)
        end
        
        Lphoton = L_photon(model, env, var)
        
        LHeat = 0.0
        if model.noneq == true
            Rate_e = Rate_volume_murca_n_e(model, core, var, model.noneq) + Rate_volume_murca_p_e(model, core, var, model.noneq)
            Rate_mu = Rate_volume_murca_n_mu(model, core, var, model.noneq) + Rate_volume_murca_p_mu(model, core, var, model.noneq)
            LHeat = var.eta_e_inf*Rate_e + var.eta_mu_inf*Rate_mu
        end

        Ce = get_Ce(model, core, var)
        Cmu = get_Cmu(model, core, var)
        Cn = get_Cn(model, core, var)
        Cp = get_Cp(model, core, var)
        
        println(ioL, "$(var.t) $(Lurca) $(LPBF) $(Lphoton) $(LHeat)")
        println(ioC, "$(var.t) $(Cn) $(Cp) $(Ce) $(Cmu)")
    end

    close(ioL)
    close(ioC)

end

function write_ini(sol, model::ModelParams)
    if isdir(model.output_dir) == false
        mkdir(model.output_dir)
    end
    
    filepath = model.output_dir * "/" * model.modelname * ".ini"
    open(filepath, "w") do f
        println(f, "# The input model parameters")
    end

    conf = ConfParse(filepath, "ini") # ini is recommended
    parse_conf!(conf)

    commit!(conf, "profile", "modelname", model.modelname)
    # starmodel
    commit!(conf, "starmodel", "eos", model.EOS)
    commit!(conf, "starmodel", "tov", model.TOV)
    commit!(conf, "starmodel", "dMoverM", model.dMoverM)
    commit!(conf, "starmodel", "del_slice", model.del_slice)
    # initial condition
    if model.noneq == true
        commit!(conf, "initial condition", "tyr0", sol[1][1])
        commit!(conf, "initial condition", "Tinf0", sol[2][1])
        commit!(conf, "initial condition", "eta_e_inf0", sol[3][1])
        commit!(conf, "initial condition", "eta_mu_inf0", sol[4][1])
    else
        commit!(conf, "initial condition", "tyr0", sol[1][1])
        commit!(conf, "initial condition", "Tinf0", sol[2][1])
        commit!(conf, "initial condition", "eta_e_inf0", 0.0)
        commit!(conf, "initial condition", "eta_mu_inf0", 0.0)
    end
    # neutron superfluidity
    commit!(conf, "neutron", "type", model.SFtype_n)
    commit!(conf, "neutron", "gap", model.gapmodel_n)
    # proton superfluidity
    commit!(conf, "proton", "type", model.SFtype_p)
    commit!(conf, "proton", "gap", model.gapmodel_p)
    # Rotochemical heating
    commit!(conf, "rotochemical", "noneq", model.noneq)
    commit!(conf, "rotochemical", "P0", model.P0)
    commit!(conf, "rotochemical", "Pnow", model.Pnow)
    commit!(conf, "rotochemical", "Pdotnow", model.Pdotnow)
    commit!(conf, "rotochemical", "Znpe", model.Znpe)
    commit!(conf, "rotochemical", "Znpmu", model.Znpmu)
    commit!(conf, "rotochemical", "Znp", model.Znp)
    commit!(conf, "rotochemical", "Wnpe", model.Wnpe)
    commit!(conf, "rotochemical", "Wnpmu", model.Wnpmu)
    # solver
    commit!(conf, "ODE", "solver", model.solver)
    commit!(conf, "ODE", "tyrf", model.tyrf)
    commit!(conf, "ODE", "reltol", model.reltol)
    commit!(conf, "ODE", "abstol", model.abstol)
    commit!(conf, "ODE", "dt", model.dt)
    # Hyper params
    commit!(conf, "Hyper params", "alpha", model.alpha)
    commit!(conf, "Hyper params", "beta", model.beta)
    # output
    commit!(conf, "output", "output_dir", model.output_dir)

    save!(conf)
end
        
end
