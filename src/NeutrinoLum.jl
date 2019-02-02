module NeutrinoLum

export L_durca_e, L_durca_mu, L_murca_n_e, L_murca_n_mu, L_murca_p_e, L_murca_p_mu

push!(LOAD_PATH, "./")

using NeutronStar
using Urca
using Simpson

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

end
