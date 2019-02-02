module ODESolvers

export cooling

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

function cooling(model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables,
                 tspan::Tuple{Float64,Float64}, reltol=1e-10, abstol=1e-10)

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
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, CVODE_Adams(), reltol, abstol)

    return sol
    
end


end
