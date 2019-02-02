module NeutrinoPBF

export Q_PBF_n, Q_PBF_p

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

using NeutronStar
using Simpson

function Q_PBF_n(T::Float64, kFn::Float64, mstn::Float64, SFtype_n::String, vn::Float64)
    """
    neutron PBF emissivity.
    Page et al., 0906.1621.
    """
        
    coeffVn = 1.0
    coeffAn = gA #g_A
    prefactor = 3.51*1e21 * mstn/mn * kFn/(mn*MeVTofminv) * (T*1e-9)^7
    
    
    if SFtype_n == "3P2m0" || SFtype_n == "3P2m2"
        a = 2*coeffAn^2
    else
        a = 0.0
        println("Q_PBF_n: $SFtype_n is not supported")
    end        
    return prefactor * a * control_func(SFtype_n, vn)
end
    

function Q_PBF_p(T::Float64, kFp::Float64, mstp::Float64, SFtype_p::String, vp::Float64)
    """
    proton PBF emissivity.
    Page et al., 0906.1621.
    """
    coeffVp = 4*sinWsq - 1.0
    coeffAp = -gA #g_A
    prefactor = 3.51*1e21 * mstp/mp * kFp/(mp*MeVTofminv) * (T*1e-9)^7
        
    if SFtype_p == "1S0"
        a = coeffVp^2 * (4.0/81.0) * (kFp/(mstp*MeVTofminv))^4 + coeffAp^2 * (kFp/(mp*MeVTofminv))^2 * (1.0 + 11.0/42.0/(mstp/mp)^2)
    else
        a = 0.0
        println("Q_PBF_p: $SFtype_p is not supported")
    end

    return prefactor * a * control_func(SFtype_p, vp)
end


function control_func(SFtype::String, v::Float64)
    """
    Fit formular for PBF control fucntion (often denoted by F or R)
    Yakovlev, Levenfish and Shibanov [astro-ph/9906456]
    """

    if SFtype == "1S0"
        x = 0.602*v^2 + 0.5942*v^4 + 0.288*v^6
        y = sqrt(0.5547 + sqrt(0.4453^2 + 0.0113*v^2))
        z = exp(-sqrt(4.0*v^2 + 2.245^2) + 2.245)
        return x*y*z

    elseif SFtype == "3P2m0"
        x = (1.204*v^2 + 3.733*v^4 + 0.3191*v^6)/(1.0 + 0.3511*v^2)
        y = (0.7591 + sqrt(0.2409^2 + 0.3145*v^2))^2
        z = exp(-sqrt(4.0*v^2 + 0.4616^2) + 0.4616)
        return x*y*z

    elseif SFtype == "3P2m2"
        return (0.4013*v^2 - 0.043*v^4 + 0.002172*v^6)/(1.0 - 0.2018*v^2 + 0.02601*v^4 - 0.001477*v^6 + 0.0000434*v^8)

    else
        println("control: $SFtype not supported")
        return 0.0
    end
end
        

end
