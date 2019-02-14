function L_photon(model::ModelParams, env::EnvelopeParams, var::StarVariables)
    return 4*pi*env.R^2 * sigmaSB * get_Teff(model, env, var)^4 * env.ephi_surface^2
end

function get_Teff(model::ModelParams, env::EnvelopeParams, var::StarVariables)
    return Teff(var.Tinf/env.ephi_surface, env.g_surface, model.dMoverM)
end

function get_Teff_inf(model::ModelParams, env::EnvelopeParams, var::StarVariables)
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
    eta = g14^2 * deltaMoverM

    Tstar = sqrt( 7.0 * Tb9 * sqrt(g14) )
    zeta = Tb9 - (Tstar*1e-3)

    # Purely iron envelope:
    # If Tb9 <x 7*10^3*sqrt(g14), zeta becomes negative.
    # In that case we just use Tstar as a crude approximation.
    if zeta > 0.0
        Teff6_Fe_4 = g14 * ( (7.0*zeta)^(2.25) + (zeta/3.0)^(1.25) )
    else
        Teff6_Fe_4 = Tstar^4
    end

    # Fully accreted envelope:
    Teff6_a_4 = g14*(18.1*Tb9)^(2.42)

    # Partially accreted envelope, interpolating the above two.
    a = (1.2 + (5.3*1e-6/eta)^(0.38)) * Tb9^(5.0/3.0)
    Teff6_4 = (a*Teff6_Fe_4 + Teff6_a_4)/(a+1.)
    return Teff6_4^(0.25) * 1e6
    # if eta  >= 1.e-30
    #     a = (1.2 + (5.3*1e-6/eta)^(0.38)) * Tb9^(5.0/3.0)
    #     Teff6_4 = (a*Teff6_Fe_4 + Teff6_a_4)/(a+1.)
    #     return Teff6_4^(0.25) * 1e6
    # else
    #     # For sufficiently small eta, 
    #     return Teff6_Fe_4^(0.25) * 1e6
    # end
        
end
