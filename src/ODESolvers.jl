module ODESolvers

export cooling, heating

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

using DifferentialEquations
using Sundials
using NeutronStar
using Setup
using HeatCapacity
using NeutrinoLum
using PhotonLum
using SuperfluidGaps
using SpinDown
#using LSODA

function cooling(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables,
                 reltol=1e-10, abstol=1e-10)
    solvers = Dict("CVODE_BDF"=>CVODE_BDF(), 
                   "CVODE_Adams"=>CVODE_Adams(),
                   "ARKODE"=>ARKODE())
    function f(u,p,t)
        var.Tinf = u
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
        return (-Lnu/C - L_photon(model, env, var)/C) * yrTosec
    end

    u0 = var.Tinf
    set_vn(model, core, var)
    set_vp(model, core, var)

    tspan = (var.t, model.tyrf)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, solvers[model.solver], reltol, abstol)

    return sol
    
end

function heating(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables,
                 reltol=1e-10, abstol=1e-10)
    
    solvers = Dict("CVODE_BDF"=>CVODE_BDF(), 
                   "CVODE_Adams"=>CVODE_Adams(),
                   "ARKODE"=>ARKODE())
    #u = [Tinf, eta_e_inf, eta_mu_inf]
    function f(du,u,p,t)
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

        du[1] = (-Lnu/C - L_photon(model, env, var)/C + var.eta_e_inf*Rate_e/C + var.eta_mu_inf*Rate_mu/C) * yrTosec
        du[2] = (-model.Znpe * Rate_e - model.Znp*Rate_mu + 2*model.Wnpe*var.Omega*var.Omega_dot) * yrTosec
        du[3] = (-model.Znp * Rate_e - model.Znpmu*Rate_mu + 2*model.Wnpmu*var.Omega*var.Omega_dot) *yrTosec
        #return [du1, du2, du3]
    end

    u0 = [var.Tinf, var.eta_e_inf, var.eta_mu_inf]
    set_vn(model, core, var)
    set_vp(model, core, var)
    set_Omega(model, var)
    set_Omega_dot(model, var)

    tspan = (var.t, model.tyrf)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, solvers[model.solver], reltol, abstol)

    return sol
    
end


end
