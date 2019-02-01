module HeatCapacity

push!(LOAD_PATH, "./")

using SpecHeat
using IntegrateData
using NeutronStar

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


end
