function get_Ce(model::ModelParams, core::StarCoreParams, var::StarVariables)
    c = spec_heat.(core.mste, core.kFe, var.Tlocal)
    #@show c
    return integrate_data(core.r_core, core.volume_elm .* c)
end

function get_Cmu(model::ModelParams, core::StarCoreParams, var::StarVariables)
    c = spec_heat.(core.mstmu, core.kFmu, var.Tlocal)
    #@show c
    return integrate_data(core.r_core, core.volume_elm .* c)
end

function get_Cn(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_n) == "normal"
        c = spec_heat.(core.mstn, core.kFn, var.Tlocal)
    else
        c = spec_heat.(core.mstn, core.kFn, var.Tlocal, var.vn, model.SFtype_n)
    end

    return integrate_data(core.r_core, core.volume_elm .* c)

end

function get_Cp(model::ModelParams, core::StarCoreParams, var::StarVariables)
    if lowercase(model.SFtype_p) == "normal"
        c = spec_heat.(core.mstp, core.kFp, var.Tlocal)
    else
        c = spec_heat.(core.mstp, core.kFp, var.Tlocal, var.vp, model.SFtype_p)
    end

    return integrate_data(core.r_core, core.volume_elm .* c)

end
