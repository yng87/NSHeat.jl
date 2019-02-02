module PhotonLum

export L_photon

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

using NeutronStar

function L_photon(model::ModelParams, env::EnvelopeParams, var::StarVariables)
    return 4*pi*env.R^2 * sigmaSB * get_Teff(model, env, var)^4 * env.ephi_surface^2
end

function get_Teff(model::ModelParams, env::EnvelopeParams, var::StarVariables)
    return Teff(var.Tinf/env.ephi_surface, env.g_surface, model.dMoverM)
end

function get_Tinf_eff(model::ModelParams, env::EnvelopeParams, var::StarVariables)
    return get_Teff(model, env, var) * env.ephi_surface
end


function Teff(Tb::Float64, g::Float64, deltaMoverM::Float64)
    """
    Non-accreted envelope
    deltaMoverM =0: Non-accreted.
    Taken from appendix A.3. of Potekhin, Chabrier, Yakovlev (1997)
    """
    
    #Tb9 = Tb / (1.0*10**9)
    Tb9 = abs(Tb * 1e-9)
    g14 = g *1e-14 #unit of g is cm s^-2
    eta = g14*deltaMoverM
    
    Tstar = sqrt( 7.0 * Tb9 * sqrt(g14) )
    zeta = Tb9 - (Tstar*1e-3)
    
    if zeta > 0.0
        Teff6_Fe_4 = g14 * ( (7.0*zeta)^(2.25) + (zeta/3.0)^(1.25) )
    else
        Teff6_Fe_4 = Tstar^4
    end
        
    Teff6_a_4 = g14*(18.1*Tb9)^(2.42)

    #return T[K]
    if eta  >= 1.e-30
        a = (1.2 + (5.3*1e-6/eta)^(0.38)) * Tb9^(5.0/3.0)
        Teff6_4 = (a*Teff6_Fe_4 + Teff6_a_4)/(a+1.)
        return Teff6_4^(0.25) * 1e6
    else
        return Teff6_Fe_4^(0.25) * 1e6
    end
        
end

end
