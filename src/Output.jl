"""
Output functions from the solution of ODE.
sol is supposed to have the form of (t, u, retcode), where
    t: arrays of t [yr]
    u: arrays of Tinf [K] or (Tinf[K], eta_e_inf[erg], eta_mu_inf[erg])
    retcode: success or not
"""

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
    filepath_R = model.output_dir * "/eta_rate.dat"

    ioL = open(filepath_L, "w")
    ioC = open(filepath_C, "w")
    ioR = open(filepath_R, "w")

    println(ioL, "# t[yr] Urca PBF Photon Heat")
    println(ioL, "# Luminosity unit = [erg/s]")
    println(ioC, "# t[yr] n p e mu")
    println(ioC, "# Heat Capacity unit = [erg/K]")
    println(ioR, "# t[yr] Znpe*Re Znp*R_mu Znp*R_e Znpmu*R_mu 2Wnpe*Omega*Omegadot 2Wnpmu*Omega*Omegadot")
    println(ioR, "# Unit = erg/s")
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
        Znpe_Re = 0.0
        Znp_Rmu = 0.0
        Znp_Re = 0.0
        Znpmu_Rmu = 0.0
        Wnpe_O = 0.0
        Wnpmu_O = 0.0
        if model.noneq == true
            Rate_e = Rate_volume_murca_n_e(model, core, var, model.noneq) + Rate_volume_murca_p_e(model, core, var, model.noneq)
            Rate_mu = Rate_volume_murca_n_mu(model, core, var, model.noneq) + Rate_volume_murca_p_mu(model, core, var, model.noneq)
            LHeat = var.eta_e_inf*Rate_e + var.eta_mu_inf*Rate_mu

            Znpe_Re = model.Znpe * Rate_e
            Znp_Rmu = model.Znp*Rate_mu
            Znp_Re = model.Znp * Rate_e
            Znpmu_Rmu = model.Znpmu*Rate_mu
            Wnpe_O = 2*model.Wnpe*var.Omega*var.Omega_dot
            Wnpmu_O = 2*model.Wnpmu*var.Omega*var.Omega_dot
        end

        Ce = get_Ce(model, core, var)
        Cmu = get_Cmu(model, core, var)
        Cn = get_Cn(model, core, var)
        Cp = get_Cp(model, core, var)
        
        println(ioL, "$(var.t) $(Lurca) $(LPBF) $(Lphoton) $(LHeat)")
        println(ioC, "$(var.t) $(Cn) $(Cp) $(Ce) $(Cmu)")
        println(ioR, "$(var.t) $(Znpe_Re) $(Znp_Rmu) $(Znp_Re) $(Znpmu_Rmu) $(Wnpe_O) $(Wnpmu_O)")
    end

    close(ioL)
    close(ioC)
    close(ioR)
    
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
    # output
    commit!(conf, "output", "output_dir", model.output_dir)
    # DM
    commit!(conf, "DM", "DM_heating", model.DM_heating)
    commit!(conf, "DM", "ann_fraction", model.ann_fraction)
    commit!(conf, "DM", "f_capture", model.f_capture)
    commit!(conf, "DM", "v_DM", model.v_DM)
    commit!(conf, "DM", "rho_DM", model.rho_DM)

    save!(conf)
end
        
