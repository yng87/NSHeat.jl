push!(LOAD_PATH, "../src/")

@show  LOAD_PATH

using Setup
using HeatCapacity
using NeutronStar
using SpecHeat
using Simpson
using SuperfluidGaps

function main()
    println(PROGRAM_FILE," start!!")

    model, core, env, var = setup("../src/sample.ini")
    
    var.Tinf = 1e8
    set_Tlocal(core, var)
    set_vn(model, core, var)
    set_vp(model, core, var)
    @show var.Tlocal
    Cs = get_Cn(model, core, var)

    c0s = spec_heat(core.mstn[1], core.kFn[1], var.Tlocal[1], var.vn[1], model.SFtype_n)

    c = spec_heat.(core.mstn, core.kFn, var.Tlocal)
    Cn = integrate_data(core.r_core, core.volume_elm .* c)

    @show Cs/Cn
    @show Cs
    @show Cn

    @show core.kFn[1]
    @show core.Tc_n[1]
    @show core.Tc_n[1] / var.Tlocal[1]
    @show c0s / c[1]
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

