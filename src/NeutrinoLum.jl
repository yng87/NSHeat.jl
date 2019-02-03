module NeutrinoLum

export L_durca_e, L_durca_mu, L_murca_n_e, L_murca_n_mu, L_murca_p_e, L_murca_p_mu, L_PBF_n, L_PBF_p,
    Rate_volume_murca_n_e, Rate_volume_murca_n_mu, Rate_volume_murca_p_e, Rate_volume_murca_p_mu

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

using NeutronStar
using Urca
using Simpson
using NeutrinoPBF
using UrcaNoneq

function L_durca_e(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_durca.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe)
    else
        q = Q_durca.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                     model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end

    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_durca_mu(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_durca.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu)
    else
        q = Q_durca.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                     model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end

    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_murca_n_e(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe)
        
    elseif lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) != "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_p, var.vp, var.Tlocal./core.Tc_p)
        
    elseif lowercase(model.SFtype_n) != "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_n, var.vn, var.Tlocal./core.Tc_n)
    else
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end
    
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_murca_n_mu(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu)
        
    elseif lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) != "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_p, var.vp, var.Tlocal./core.Tc_p)
        
    elseif lowercase(model.SFtype_n) != "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_n, var.vn, var.Tlocal./core.Tc_n)
    else
        q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_murca_p_e(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe)
        
    elseif lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) != "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_p, var.vp)
        
    elseif lowercase(model.SFtype_n) != "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_n, var.vn)
    else
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                       model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end
    
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)

end

function L_murca_p_mu(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu)
        
    elseif lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) != "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_p, var.vp)
        
    elseif lowercase(model.SFtype_n) != "normal" && lowercase(model.SFtype_p) == "normal"
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_n, var.vn)
    else
        q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                       model.SFtype_n, model.SFtype_p, var.vn, var.vp)
    end
    
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_PBF_n(model::ModelParams, core::StarCoreParams, var::StarVariables)
    q = Q_PBF_n.(var.Tlocal, core.kFn, core.mstn, model.SFtype_n, var.vn)
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

function L_PBF_p(model::ModelParams, core::StarCoreParams, var::StarVariables)
    q = Q_PBF_p.(var.Tlocal, core.kFp, core.mstp, model.SFtype_p, var.vp)
    return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
end

""" 
Non-equilibrium
"""
function L_murca_n_e(model::ModelParams, core::StarCoreParams, var::StarVariables,
                     eq::Bool)
    if eq == false
        xi = var.eta_e_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe, xi)
        else
            q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                           model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                           xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
    else
        return L_murca_n_e(model, core, var)
    end
    
end

function L_murca_n_mu(model::ModelParams, core::StarCoreParams, var::StarVariables,
                      eq::Bool)
    if eq == false
        xi = var.eta_mu_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu, xi)
        else
            q = Q_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                           model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                           xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
    else
        return L_murca_n_mu(model, core, var)
    end
end


function Rate_volume_murca_n_e(model::ModelParams, core::StarCoreParams, var::StarVariables,
                               eq::Bool)
     if eq == false
        xi = var.eta_e_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Rate_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe, xi)
        else
            q = Rate_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                              model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                              xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi)
    else
         return 0.0
     end
end

function Rate_volume_murca_n_mu(model::ModelParams, core::StarCoreParams, var::StarVariables,
                          eq::Bool)
    if eq == false
        xi = var.eta_mu_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Rate_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu, xi)
        else
            q = Rate_murca_n.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                              model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                              xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi)
    else
         return 0.0
     end
end

function L_murca_p_e(model::ModelParams, core::StarCoreParams, var::StarVariables,
                     eq::Bool)
    if eq == false
        xi = var.eta_e_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe, xi)
        else
            q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                           model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                           xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
    else
        return L_murca_p_e(model, core, var)
    end
    
end

function L_murca_p_mu(model::ModelParams, core::StarCoreParams, var::StarVariables,
                      eq::Bool)
    if eq == false
        xi = var.eta_mu_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu, xi)
        else
            q = Q_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                           model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                           xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi.^2)
    else
        return L_murca_p_mu(model, core, var)
    end
end


function Rate_volume_murca_p_e(model::ModelParams, core::StarCoreParams, var::StarVariables,
                               eq::Bool)
     if eq == false
        xi = var.eta_e_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Rate_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe, xi)
        else
            q = Rate_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mste, core.kFn, core.kFp, core.kFe,
                              model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                              xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi)
    else
         return 0.0
     end
end

function Rate_volume_murca_p_mu(model::ModelParams, core::StarCoreParams, var::StarVariables,
                                eq::Bool)
    if eq == false
        xi = var.eta_mu_inf / (kB*var.Tinf)
        if lowercase(model.SFtype_n) == "normal" && lowercase(model.SFtype_p) == "normal"
            q = Rate_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu, xi)
        else
            q = Rate_murca_p.(var.Tlocal, core.mstn, core.mstp, core.mstmu, core.kFn, core.kFp, core.kFmu,
                              model.SFtype_n, model.SFtype_p, var.vn, var.vp,
                              xi)
        end
        return integrate_data(core.r_core, core.volume_elm .* q .* core.ephi)
    else
         return 0.0
     end
end

end
