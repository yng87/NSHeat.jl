push!(LOAD_PATH, "../src/")

@show  LOAD_PATH

using Setup
using HeatCapacity
using NeutronStar
using SpecHeat
using Simpson
using SuperfluidGaps
using NeutrinoLum
using DifferentialEquations

function main()
    println(PROGRAM_FILE," start!!")

    model, core, env, var = setup("../src/sample.ini")
    
    var.Tinf = 1e8
    set_Tlocal(core, var)
    set_vn(model, core, var)
    set_vp(model, core, var)
    @show var.Tlocal
    @show var.vn
    @show var.vp
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

    println("##########")
    @show L_durca_e(model, core, var)
    @show L_durca_mu(model, core, var)
    @show L_murca_n_e(model, core, var)
    @show L_murca_n_mu(model, core, var)
    @show L_murca_p_e(model, core, var)
    @show L_murca_p_mu(model, core, var)

    function f(u,p,t)
        var.Tinf = u
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)

        C = get_Ce(model, core, var) + get_Cmu(model, core, var) + get_Cn(model, core, var) + get_Cp(model, core, var)
        Lnu = L_murca_n_e(model, core, var) + L_murca_n_mu(model, core, var) + L_murca_p_e(model, core, var) + L_murca_p_mu(model, core, var)
        return -Lnu/C
    end
    var.Tinf = 1e10
    u0 = var.Tinf
    @show u0
    tspan = (0.1, 1e5)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, alg_hints=[:stiff], reltol=1e-8, abstol=1e-8)

    @show sol
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

