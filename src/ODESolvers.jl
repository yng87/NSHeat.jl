"""
Functions whicn integrate all the parts like heat capacity or luminosity, 
and solve the 1-D (spherically symmetric and isothermal) temperature development.

If you add another functions which affects cooling/heating, do not forget to modity these function properly.

Among the varous solvers, CVODE_BDF is the best.
But for some parameter choice, the ODE becomes really stiff and even CVODE_BDF works not so well.
In such a case, you must tune the tolerance you impose on the ODE solver.
"""

function cooling(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)

    solvers = Dict("CVODE_BDF"=>CVODE_BDF(), 
                   "CVODE_Adams"=>CVODE_Adams(),
                   "ARKODE"=>ARKODE(),
                   "Rosenbrock23"=>Rosenbrock23(autodiff=false),
                   "TRBDF2"=>TRBDF2(autodiff=false),
                   "Rodas5"=>Rodas5(autodiff=false),
                   "Rodas4P"=>Rodas4P(autodiff=false),
                   "Kvaerno5"=>Kvaerno5(autodiff=false),
                   "KenCarp4"=>KenCarp4(autodiff=false),
                   "radau"=>radau())
    
    function f(u,p,t)
        model, core, env, var = p
        var.t = exp(t) #yr
        var.Tinf = exp(u)
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)
        #Heat capacity
        C = get_Ce(model, core, var) + get_Cmu(model, core, var) + get_Cn(model, core, var) + get_Cp(model, core, var)
        #Neutrino luminosity
        Lnu = L_murca_n_e(model, core, var) + L_murca_n_mu(model, core, var) + L_murca_p_e(model, core, var) + L_murca_p_mu(model, core, var)
        if lowercase(model.SFtype_n) != "normal"
            Lnu += L_PBF_n(model, core, var)
        end
        if lowercase(model.SFtype_p) != "normal"
            Lnu += L_PBF_p(model, core, var)
        end
        #Do not forget yrTosec!
        du = (-Lnu/C - L_photon(model, env, var)/C) * yrTosec * var.t/var.Tinf
        return du
    end

    u0 = log(var.Tinf)
    set_vn(model, core, var)
    set_vp(model, core, var)

    tspan = (log(var.t), log(model.tyrf))
    p = (model, core, env, var)
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, solvers[model.solver], reltol=model.reltol, abstol=model.abstol, saveat=model.dt)

    return exp.(sol.t), exp.(sol[1,:]), sol.retcode
    
end

function heating(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    
    solvers = Dict("CVODE_BDF"=>CVODE_BDF(linear_solver=:GMRES, max_convergence_failures=1000), 
                   "CVODE_Adams"=>CVODE_Adams(),
                   "ARKODE"=>ARKODE(linear_solver=:GMRES),
                   "Rosenbrock23"=>Rosenbrock23(autodiff=false),
                   "TRBDF2"=>TRBDF2(autodiff=false),
                   "Rodas5"=>Rodas5(autodiff=false),
                   "Rodas4P"=>Rodas4P(autodiff=false),
                   "Kvaerno5"=>Kvaerno5(autodiff=false),
                   "KenCarp4"=>KenCarp4(autodiff=false),
                   "radau"=>radau(),
                   "lsoda"=>lsoda())

    #u = [Tinf, eta_e_inf, eta_mu_inf]
    function f(du,u,p,t)
        # u = ln([T, eta_e_inf, eta_mu_inf])
        # t = ln(t/yr)
        model, core, env, var = p
        var.t = exp(t) #yr
        var.Tinf = exp(u[1])
        var.eta_e_inf = exp(u[2]) #erg
        var.eta_mu_inf = exp(u[3]) #erg
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)
        set_Omega(model, var)
        set_Omega_dot(model, var)
        #Heat capacity
        C = get_Ce(model, core, var) + get_Cmu(model, core, var) + get_Cn(model, core, var) + get_Cp(model, core, var)
        #Neutrino luminosity
        Lnu = L_murca_n_e(model, core, var, model.noneq) + L_murca_n_mu(model, core, var, model.noneq) + L_murca_p_e(model, core, var, model.noneq) + L_murca_p_mu(model, core, var, model.noneq)
        if lowercase(model.SFtype_n) != "normal"
            Lnu += L_PBF_n(model, core, var)
        end
        if lowercase(model.SFtype_p) != "normal"
            Lnu += L_PBF_p(model, core, var)
        end
        #Do not forget yrTosec!
        Rate_e = Rate_volume_murca_n_e(model, core, var, model.noneq) + Rate_volume_murca_p_e(model, core, var, model.noneq)
        Rate_mu = Rate_volume_murca_n_mu(model, core, var, model.noneq) + Rate_volume_murca_p_mu(model, core, var, model.noneq)

        du[1] = (-Lnu/C - L_photon(model, env, var)/C + var.eta_e_inf*Rate_e/C + var.eta_mu_inf*Rate_mu/C) * yrTosec * var.t/var.Tinf
        du[2] = (-model.Znpe * Rate_e - model.Znp*Rate_mu + 2*model.Wnpe*var.Omega*var.Omega_dot) * yrTosec * var.t/var.eta_e_inf
        du[3] = (-model.Znp * Rate_e - model.Znpmu*Rate_mu + 2*model.Wnpmu*var.Omega*var.Omega_dot) *yrTosec * var.t/var.eta_mu_inf
        #return du
    end

    u0 = log.([var.Tinf, var.eta_e_inf, var.eta_mu_inf])
    set_vn(model, core, var)
    set_vp(model, core, var)
    set_Omega(model, var)
    set_Omega_dot(model, var)

    tspan = (log(var.t), log(model.tyrf))
    p = (model, core, env, var)
    prob = ODEProblem(f, u0, tspan, p)

    sol = solve(prob, solvers[model.solver], abstol=model.abstol, reltol=model.reltol, saveat=model.dt)

    return exp.(sol.t), exp.(sol[1,:]), exp.(sol[2,:]), exp.(sol[3,:]), sol.retcode
    
end

function heating_lin(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    
    solvers = Dict("CVODE_BDF"=>CVODE_BDF(linear_solver=:GMRES, max_convergence_failures=1000), 
                   "CVODE_Adams"=>CVODE_Adams(),
                   "ARKODE"=>ARKODE(linear_solver=:GMRES),
                   "Rosenbrock23"=>Rosenbrock23(autodiff=false),
                   "TRBDF2"=>TRBDF2(autodiff=false),
                   "Rodas5"=>Rodas5(autodiff=false),
                   "Rodas4P"=>Rodas4P(autodiff=false),
                   "Kvaerno5"=>Kvaerno5(autodiff=false),
                   "KenCarp4"=>KenCarp4(autodiff=false),
                   "radau"=>radau(),
                   "lsoda"=>lsoda())

    #u = [Tinf, eta_e_inf, eta_mu_inf]
    function f(du,u,p,t)
        # t = ln(t/yr)
        model, core, env, var = p
        var.t = t #yr
        var.Tinf = u[1]
        var.eta_e_inf = u[2] #erg
        var.eta_mu_inf = u[3] #erg
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)
        set_Omega(model, var)
        set_Omega_dot(model, var)
        #Heat capacity
        C = get_Ce(model, core, var) + get_Cmu(model, core, var) + get_Cn(model, core, var) + get_Cp(model, core, var)
        #Neutrino luminosity
        Lnu = L_murca_n_e(model, core, var, model.noneq) + L_murca_n_mu(model, core, var, model.noneq) + L_murca_p_e(model, core, var, model.noneq) + L_murca_p_mu(model, core, var, model.noneq)
        if lowercase(model.SFtype_n) != "normal"
            Lnu += L_PBF_n(model, core, var)
        end
        if lowercase(model.SFtype_p) != "normal"
            Lnu += L_PBF_p(model, core, var)
        end
        #Do not forget yrTosec!
        Rate_e = Rate_volume_murca_n_e(model, core, var, model.noneq) + Rate_volume_murca_p_e(model, core, var, model.noneq)
        Rate_mu = Rate_volume_murca_n_mu(model, core, var, model.noneq) + Rate_volume_murca_p_mu(model, core, var, model.noneq)
        @show t, u
        du[1] = (-Lnu/C - L_photon(model, env, var)/C + var.eta_e_inf*Rate_e/C + var.eta_mu_inf*Rate_mu/C) * yrTosec
        du[2] = (-model.Znpe * Rate_e - model.Znp*Rate_mu + 2*model.Wnpe*var.Omega*var.Omega_dot) * yrTosec 
        du[3] = (-model.Znp * Rate_e - model.Znpmu*Rate_mu + 2*model.Wnpmu*var.Omega*var.Omega_dot) *yrTosec 
        #return du
    end

    u0 = [var.Tinf, var.eta_e_inf, var.eta_mu_inf]
    set_vn(model, core, var)
    set_vp(model, core, var)
    set_Omega(model, var)
    set_Omega_dot(model, var)

    tspan = (var.t, model.tyrf)
    p = (model, core, env, var)
    prob = ODEProblem(f, u0, tspan, p)
    saveat = exp.(log(var.t):model.dt:log(model.tyrf))
    sol = solve(prob, solvers[model.solver], abstol=model.abstol, reltol=model.reltol, saveat=saveat)

    return sol.t, sol[1,:], sol[2,:], sol[3,:], sol.retcode
    
end


