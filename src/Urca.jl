module Urca

export Q_murca_n, Q_murca_p

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

function Q_durca()
    return
end

function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64)
    #e [erg/s/cm^3]
    T9 = T * 1e-9
    k0 = 1.68 # corresponds to n0 = 0.16  fm^-3
    alpha = 1.76 - 0.63*(k0/kn)^2
    beta = 0.68
    vFl = pFl/mstl # Fermi velocity for lepton
    return 8.05e21 * (mstn/mn)^3 * (mstp/mp) * (kFp/k0) * T9^8 *alpha*beta * vFl
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64)
    
    return elseif(3*kFp+kFl-kFn<0,
                  0,
                  Q_murca_n(T, mstn, mstp, mstl, kFn, kFp) * (mstp/mstn)^3 * (kFl+3*kFp-kFn)^2/(8*kFl*kFp)
                  )
end

function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64,
                   SFtype::String, v::Float64, t::Float64)
    if SFtype == "1S0"
        # assume proton superfluidity
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp) * Rn_SFp(v)
    elseif SFtype == "3P2m0"
        # argument ot R is t=T/Tc, not v
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp) * Rn_SFn(t)
    else
        println("Q_murca_n: SF type not supported")
        return 0.0
    end
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype::String, v::Float64)
    if SFtype == "1S0"
        # assume proton superfluidity
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl) * Rp_SFp(v)
    elseif SFtype == "3P2m0"
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl) * Rn_SFn(v)
    else
        println("Q_murca_n: SF type not supported")
        return 0.0
    end
end

function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64)
    # only for proton 1S0 and neutron 3P2m0
    return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp) * Rn_SFnp(vn, vp)
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64)
    # only for proton 1S0 and neutron 3P2m0
    return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl) * Rp_SFnp(vn, vp)
    
"""
Superfluid reduction factor for modified Urca process.
Numerical integration is done with Gauss-Laguerre quadrature and trapezotal rule.

Reduction factor is named I{n,p}_SF{n,p}{n,p}_fit, meaning that
- In (Ip): reduction factor for neutron (proton) branch.
- SFn (SFp): only neutron (proton) superfluidity is turned on 
- SFnp: both proton and neutron are superfluid.

The fitting formulas are taken from 
Yakovlev et al, Physics reports Volume 354, Issues 1–2, November 2001, Pages 1-155
Gusakov Gusakov, A&A 389, 702-715 (2002)
"""

"""
Only proton superfluidity
"""

function Rn_SFp(vp::Float64)
    a = 0.1477 + sqrt(0.8523^2 + (0.1175*vp)^2)
    b = 0.1477 + sqrt(0.8523^2 + (0.1297*vp)^2)
    return (a^7.5 + b^5.5)/2.0 * exp(3.437 - sqrt(3.437^2 + vp^2))
end

function Rp_SFp(vp::Float64)
    return (0.2414 + sqrt(0.7586^2 + (0.1318*vp)^2))^7 * exp(5.339 - sqrt(5.339^2 + 4*vp^2))
end


"""
Only neutron superfluidity
"""

vB(t) = sqrt(1-t)*(0.7893 + 1.188/t)
function Rn_SFn(tn::Float64)
    # coming from NSCool
    # I changed in prefactor 39.1 to 3.91 just to obtain better fitting
    vn = vB(tn)
    c = 0.1612 + sqrt(0.8388^2 + (0.1117*vn)^2)
    d = 0.1612 + sqrt(0.8388^2 + (0.1274*vn)^2)
    prefactor = 3.91*tn*exp(-1.188/tn)
    return prefactor * (c^7+d^5)/2.0 * exp(2.398 - sqrt(2.398^2 + vn^2))
end

function Rp_SFn(vn::Float64)
     a = 0.1612 + sqrt(0.8388^2 + (0.1117*vn)^2)
     b = 0.1612 + sqrt(0.8388^2 + (0.1274*vn)^2)
     return (a^7+b^5)/2.0 * exp(2.398 - sqrt(2.398^2 + vn^2))
end


"""
Both proton and neutron superfluidity
"""
pn = Array{Float64,1}[]
open("../fit_table/Murca_reduction_nbr_AB_fit_coeffs.dat", "r") do f
    for line in eachline(f)
        if line[1]!='#'
            push!(pn, (parse.(Float64, split(line, ","))))
        end
    end
end

pp = Array{Float64,1}[]
open("../fit_table/Murca_reduction_pbr_AB_fit_coeffs.dat", "r") do f
    for line in eachline(f)
        if line[1]!='#'
            push!(pp, (parse.(Float64, split(line, ","))))
        end
    end
end

function Rn_SFnp(v1::Float64, v2::Float64)
    ps = pn
    # v1 = vn, v2 = vp
    v = sqrt(v1^2 + v2^2) 
    phi = asin(v1/v)
    
    # Region 4
    if v <= 5.0
        p = ps[4]
        
        A = p[1]*v1^2 + p[2]*v2^2 + p[3]*v1^2*v2^2 + p[4]*v1^6
        B = sqrt(1 + p[6]*v1^2 + p[7]*v2^2 + p[5]*v1^8 + p[8]*v1^6)
        C = 1 + p[9]*v2^4
        
        return C * exp(-A/B)
    
    # Region 1
    elseif v1 > v2
        p = ps[1]
        y = sin(phi+p[15])^2
        t = cos(phi)^2
        z = cos(phi+p[14])^2
        
        A = p[1] + p[2]*t^2 + p[4]/(1+p[3]*t) - p[5]*t
        B = p[6] + p[7]*t^2 + p[9]*t/(1+p[8]*t^2) - p[10]*t
        C = p[11] - p[12]/y/(1+p[13]*z^2)^3
        
        return exp(-A*v^2/(1+B*v^2)^C)
        
    # Region 2
    elseif v2>=v1 && v1>=v2/3.0
        p = ps[2]
        z = cos(phi)^2
        t = cos(phi+p[8])^2
        q  = sin(p[9]*phi+p[10])^2
        y = sin(phi+3*pi/4.0)^2
        
        A = p[1] + p[2]*z^2 + p[4]/(1+p[3]*z) + p[5]*z
        B = 0.035
        C = p[6] + p[7]*z^2 + p[9]/(1+p[8]*z) + p[10]*z
        
        return exp(-A*v^2/(1+B*v^2)^C)
        
    # Region 3
    elseif v1 < v2/3.0
        p = ps[3]
        A = p[12]*(1 + p[15]*phi^2 + p[16]*phi^3 + p[17]*phi^4)/(1 + p[13]*phi^2 + p[14]*phi^3)
        B = p[11] + p[7]/(1 + p[8]*phi^2 + p[9]*phi^3 + p[10]*phi^4)
        C = p[1]*(1 + p[4]*phi^2 + p[5]*phi^3 + p[6]*phi^4)/(1 + p[2]*phi^2 + p[3]*phi^3)
        
        return exp(-A*v^2/(1+B*v^2)^C)
    end
end

function Rp_SFnp(v1::Float64, v2::Float64)
    ps = pp
    # v1 = vn, v2 = vp
    v = sqrt(v1^2 + v2^2) 
    phi = asin(v1/v)
    
    # Region 4
    if v <= 5.0
        p = ps[4]
        
        A = p[1]*v2^2 + p[2]*v1^2 + p[3]*v1^2*v2^2 + p[4]*v2^6 + p[5]*v1^6
        B = 1 + p[6]*v2^2 + p[7]*v1^2 + p[8]*v2^4
        C = 1 + p[9]*v1^4
        
        return C * exp(-A/B)
    
    # Region 1
    elseif v1 > v2
        p = ps[1]
        y = sin(phi+p[15])^2
        t = cos(phi)^2
        
        A = p[1] + p[2]*phi + p[4]*phi/(1+p[3]*t*phi)^2 + p[5]*t*phi^2
        B = p[6] + p[7]*phi + p[9]*phi/(1+p[8]*t*phi+p[10]*y*t)^2 + p[11]*phi^2
        C = p[12] - p[13]*t + p[16]/(1+p[14]*t^2)^2 + p[17]*t*phi
        
        return exp(-A*v^2/(1+B*v^2)^C)
        
    # Region 2
    elseif v2>=v1 && v1>=v2/3.0
        p = ps[2]
        z = cos(phi)^2
        
        A = p[1] + p[2]*z^2 + p[4]/(1+p[3]*z) + p[5]*phi
        B = p[6] + p[7]*z^2 + p[8]/(1+p[9]*z) + p[10]*phi
        C = p[11] - p[12]*phi^2 + p[13]*phi^2/(1+p[14]*z) + p[15]*phi
        
        return exp(-A*v^2/(1+B*v^2)^C)
        
    # Region 3
    elseif v1 < v2/3.0
        p = ps[3]
        A = p[12]*(1 + p[15]*phi^2 + p[16]*phi^3 + p[17]*phi^4)/(1 + p[13]*phi^2 + p[14]*phi^3)
        B = p[11] + p[7]/(1 + p[8]*phi^2 + p[9]*phi^3 + p[10]*phi^4)
        C = p[1]*(1 + p[4]*phi^2 + p[5]*phi^3 + p[6]*phi^4)/(1 + p[2]*phi^2 + p[3]*phi^3)
        
        return exp(-A*v^2/(1+B*v^2)^C)
    end
end

end