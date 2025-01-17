"""
Functions for superfluid gap.

A : 1S0
B : 3P2(mJ=0)
C : 3P2(mJ=2)

"""

"""
Gap functions
Yakovlev et. al, Phys.Rept. 354 (2001) 1
"""
vA(t) = sqrt(1-t)*(1.456 - 0.157/sqrt(t) + 1.764/t)
vB(t) = sqrt(1-t)*(0.7893 + 1.188/t)
vC(t) = sqrt(1-t^4)/t*(2.03 - 0.4903*t^4 + 0.1727*t^8)

function fit_frac(kF::Float64, a::Float64, k1::Float64, k2::Float64, k3::Float64, k4::Float64)
    return ifelse( kF>= k1 && kF < k3, a * (kF-k1)^2/((kF-k1)^2 + k2^2) * (kF-k3)^2/((kF-k3)^2 + k4^2) *1e10, 0.0)
end

function fit_gauss(kF::Float64, kfmax::Float64, delkf::Float64, tcmax::Float64)
    return tcmax * exp(-(kF-kfmax)^2/delkf^2)
end

function fit_mod_gauss(kF::Float64, kfmax::Float64, delkf::Float64, tcmax::Float64, r4::Float64)
    return tcmax * exp(-(kF-kfmax)^2/delkf^2 - r4 * (kF-kfmax)^4/delkf^4)
end

gap_params_p_S = Dict("AO"=>(4.98, 0.038, 0.648, 1.05, 1.41),
                      "T73"=>(0.845, 0.208, 0.276, 0.844, 0.245),
                      "CCDK"=>(26.6933, 0.0206134, 1.67406, 1.27664, 1.26246),
                      "AO_mod_gauss"=>(0.497, 0.293, 2.33e9, 0.0253),
                      "CCDK_mod_gauss"=>(0.66, 0.457, 6.50e9, 0.636),
                      "BS_2BF_mod_gauss"=>(0.640, 0.558, 7.97e9, 1.54),
                      "BS_2BF_3BF_mod_gauss"=>(0.477, 0.338, 5.67e9, 0.434),
                      "EEHO_mod_gauss"=>(0.624, 0.422, 5.73e9, 0.771)
                      )
        
gap_params_n_Pm0 = Dict("a"=>(1.8, 0.5, 1e9),
                        "b"=>(2.0, 0.5, 3e9),
                        "c"=>(2.5, 0.7, 1e10),
                        "a2"=>(2.3, 0.9, 5.5e8),
                        "H"=>(4.8*MeVToerg/1.188/kB*1e-10, 1.07, sqrt(1.8), 3.2, sqrt(2.0))
                        )

function Tc_p_S(model_name::String, kF::Float64)
    # We can extend model space, e.g. including modified gaussian
    if occursin("mod_gauss", model_name)
        return fit_mod_gauss(kF, gap_params_p_S[model_name]...)
    else
        return fit_frac(kF, gap_params_p_S[model_name]...)
    end
end

function Tc_n_Pm0(model::ModelParams, kF::Float64)
    model_name = model.gapmodel_n
    if model_name == "H"
        return fit_frac(kF, gap_params_n_Pm0[model_name]...)
    elseif model_name == "mod_gauss_custom"
        return fit_mod_gauss(kF, model.kfmax_n, model.delkf_n, model.tcmax_n, model.r4_n)
    else
        return fit_gauss(kF, gap_params_n_Pm0[model_name]...)
    end

end

function set_Tc_p(model::ModelParams, kF::Array{Float64,1})
    if model.SFtype_p == "1S0"
        Tc_p = Tc_p_S.(model.gapmodel_p, kF)
    elseif lowercase(model.SFtype_p) == "normal"
        Tc_p = zeros(length(kF))
    else
        Tc_p = zeros(length(kF))
        println("The proton gap type specificatoin is wrong, Tc is set to be zero")
    end
    
    return Tc_p
end

function set_Tc_n(model::ModelParams, kF::Array{Float64,1})
    if model.SFtype_n == "3P2m0"
        #Tc_n = Tc_n_Pm0.(model, kF)
        Tc_n = map(x->Tc_n_Pm0(model, x), kF)
    elseif model.SFtype_n == "3P2m2"
        Tc_n = zeros(length(kF))
        println("Neutron 3P2(mJ=2) gap is not implemented yet!, Tc is set to be zero")
    elseif lowercase(model.SFtype_n) == "normal"
        Tc_n = zeros(length(kF))
    else
        Tc_n = zeros(length(kF))
        println("The neutron gap type specificatoin is wrong, Tc is set to be zero")
    end

    return Tc_n
end


# The following two functions should be firther optimized
function set_vn(model::ModelParams, core::StarCoreParams, var::StarVariables)
    sf_idx = var.Tlocal .< core.Tc_n
    t = var.Tlocal[sf_idx] ./ core.Tc_n[sf_idx]

    # To avoid negative t. Just tempral, should make callback function.
    neg_t = .!(t .> 0)
    t[neg_t] = zeros(length(neg_t))[neg_t] .+ 1e-10
        
    if model.SFtype_n == "3P2m0"
        var.vn[sf_idx] = vB.(t)
        var.vn[.!sf_idx] = zeros(length(var.vn))[.!sf_idx]
    elseif model.SFtype_n == "3P2m2"
        println("set_vn: neutron 3P2m2 is not supported")
    elseif lowercase(model.SFtype_n) == "normal"
        var.vn = zeros(length(var.vn))
    else
        println("set_vn: wrong neutron superfluid type")
    end
    
    return nothing
end

function set_vp(model::ModelParams, core::StarCoreParams, var::StarVariables)
    sf_idx = var.Tlocal .< core.Tc_p
    t = var.Tlocal[sf_idx] ./ core.Tc_p[sf_idx]

    # To avoid negative t. Just tempral, should make callback function.
    neg_t = .!(t .> 0)
    t[neg_t] = zeros(length(neg_t))[neg_t] .+ 1e-10

    if model.SFtype_p == "1S0"
        var.vp[sf_idx] = vA.(t)
        var.vp[.!sf_idx] = zeros(length(var.vp))[.!sf_idx]
    elseif lowercase(model.SFtype_p) == "normal"
        var.vp = zeros(length(var.vp))
    else
        println("set_vn: wrong neutron superfluid type")
    end
    return nothing
end
