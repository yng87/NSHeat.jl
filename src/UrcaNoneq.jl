module UrcaNoneq

export Q_durca, Q_murca_n, Q_murca_p, Rate_durca, Rate_murca_n, Rate_murca_p

push!(LOAD_PATH, "/")
include("./PhysicalConstants.jl")

#using Urca
import Urca:Q_durca, Q_murca_n, Q_murca_p
using Dierckx
using DelimitedFiles
using MurcaNoneqNumerical

"""
Non-superfluid
"""

function Q_durca(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                 xi::Float64)
    return Q_durca(T, mstn, mstp, mstl, kFn, kFp, kFl) * FD(xi)
end

function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   xi::Float64)
    return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl) * FM(xi)
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   xi::Float64)
    return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl) * FM(xi)
end

function Rate_durca(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                    xi::Float64)
    return Q_durca(T, mstn, mstp, mstl, kFn, kFp, kFl)/(kB*T) * HD(xi)
end

function Rate_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      xi::Float64)
    return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl)/(kB*T) * HM(xi)
end

function Rate_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      xi::Float64)
    return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl)/(kB*T) * HM(xi)
end

"""
Superfluid
"""
function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                   xi::Float64)
    threshold = 1.0
    vth = 3*vn + vp
    n = 5
    xi_th = 10.0
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl,
                         SFtype_n, SFtype_p, vn, vp)
    elseif abs(xi) < vth || abs(xi) < xi_th
        # superfluid and non-equilibrium, but below threshold
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Iemis_n_SFnp(vn, vp, xi, n) / FM(xi)
    else
        # superfluid, non-equilibrium, and above threshold
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Remis_murca_n(vn, vp, xi)
    end
    
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                   xi::Float64)
    threshold = 1.0
    vth = vn + 3*vp
    n = 5
    xi_th = 10.0
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl,
                         SFtype_n, SFtype_p, vn, vp)
    elseif abs(xi) < vth || abs(xi) < xi_th
        # superfluid and non-equilibrium, but below threshold
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Iemis_p_SFnp(vn, vp, xi, n) / FM(xi)
    else
        # superfluid, non-equilibrium, and above threshold
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Remis_murca_p(vn, vp, xi)
    end
    
end

function Rate_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                      xi::Float64)

    threshold = 1.0
    vth = 3*vn + vp
    n = 5
    xi_th = 10.0
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Rate_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return 0.0
    elseif abs(xi) < vth || abs(xi) < xi_th
        # superfluid and non-equilibrium, but below threshold
        return Rate_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Irate_n_SFnp(vn, vp, xi, n) / HM(xi)
    else
        # superfluid, non-equilibrium, and above threshold
        return Rate_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Rrate_murca_n(vn, vp, xi)
    end

end

function Rate_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                      xi::Float64)
    threshold = 1.0
    vth = vn + 3*vp
    n = 5
    xi_th = 10.0
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Rate_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return 0.0
    elseif abs(xi) < vth || abs(xi) < xi_th
        # superfluid and non-equilibrium, but below threshold
        return Rate_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Irate_p_SFnp(vn, vp, xi, n) / HM(xi)
    else
        # superfluid, non-equilibrium, and above threshold
        return Rate_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Rrate_murca_p(vn, vp, xi)
    end
end

"""
Phase space / reduction factors
"""

Remis_murca_n_table = readdlm("../number_table/Remis_murca_n.dat", Float64, comments=true)
Remis_murca_n_xarr_table = readdlm("../number_table/Remis_murca_n_xarray.dat", Float64, comments=true)
Remis_murca_n_yarr_table = readdlm("../number_table/Remis_murca_n_yarray.dat", Float64, comments=true)

Remis_murca_p_table = readdlm("../number_table/Remis_murca_p.dat", Float64, comments=true)
Remis_murca_p_xarr_table = readdlm("../number_table/Remis_murca_p_xarray.dat", Float64, comments=true)
Remis_murca_p_yarr_table = readdlm("../number_table/Remis_murca_p_yarray.dat", Float64, comments=true)

Rrate_murca_n_table = readdlm("../number_table/Rrate_murca_n.dat", Float64, comments=true)
Rrate_murca_n_xarr_table = readdlm("../number_table/Rrate_murca_n_xarray.dat", Float64, comments=true)
Rrate_murca_n_yarr_table = readdlm("../number_table/Rrate_murca_n_yarray.dat", Float64, comments=true)

Rrate_murca_p_table = readdlm("../number_table/Rrate_murca_p.dat", Float64, comments=true)
Rrate_murca_p_xarr_table = readdlm("../number_table/Rrate_murca_p_xarray.dat", Float64, comments=true)
Rrate_murca_p_yarr_table = readdlm("../number_table/Rrate_murca_p_yarray.dat", Float64, comments=true)

Remis_murca_n_spl = Spline2D(Remis_murca_n_xarr_table[1,:], Remis_murca_n_yarr_table[1,:], Remis_murca_n_table)
Remis_murca_p_spl = Spline2D(Remis_murca_p_xarr_table[1,:], Remis_murca_p_yarr_table[1,:], Remis_murca_p_table)
Rrate_murca_n_spl = Spline2D(Rrate_murca_n_xarr_table[1,:], Rrate_murca_n_yarr_table[1,:], Rrate_murca_n_table)
Rrate_murca_p_spl = Spline2D(Rrate_murca_p_xarr_table[1,:], Rrate_murca_p_yarr_table[1,:], Rrate_murca_p_table)

function Remis_murca_n(vn::Float64, vp::Float64, xi::Float64)
    x = vn/xi
    y = vp/xi
    return Remis_murca_n_spl(x, y)
end

function Remis_murca_p(vn::Float64, vp::Float64, xi::Float64)
    x = vn/xi
    y = vp/xi
    return Remis_murca_p_spl(x, y)
end

function Rrate_murca_n(vn::Float64, vp::Float64, xi::Float64)
    x = vn/xi
    y = vp/xi
    return Rrate_murca_n_spl(x, y)
end

function Rrate_murca_p(vn::Float64, vp::Float64, xi::Float64)
    x = vn/xi
    y = vp/xi
    return Rrate_murca_p_spl(x, y)
end

function FM(xi::Float64)
    return 1.0 + (22020.0*xi^2)/(11513.0*pi^2) + (5670.0*xi^4)/(11513.0*pi^4) + (420.0*xi^6)/(11513.0*pi^6) + (9.0*xi^8)/(11513.0*pi^8)
end

function HM(xi::Float64)
    #=
    Fernandez and Reisenegger (2005), apj 625, 291-306.
    =#
    return (14680.0*xi)/(11513.0*pi^2) + (7560.0*xi^3)/(11513.0*pi^4) + (840*xi^5)/(11513.0*pi^6) + (24.0*xi^7)/(11513.0*pi^8)
end

end
