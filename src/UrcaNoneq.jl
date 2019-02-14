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

threshold = 1.0

stepfunc(xi::Float64, vth::Float64, alpha::Float64) = 1.0 / (exp(alpha*(-xi+vth)) + 1.0)

function Q_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                   xi::Float64, alpha=10.0, beta=100.0)
    vth = 3*vn + vp

    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl,
                         SFtype_n, SFtype_p, vn, vp)
    elseif abs(xi) < vth/beta
        return 0.0
    else
        return Q_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Remis_murca_n(vn, vp, xi) * stepfunc(xi, vth, alpha)
    end
    
end

function Q_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                   SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                   xi::Float64, alpha=10.0, beta=100.0)
    vth = vn + 3*vp
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl,
                         SFtype_n, SFtype_p, vn, vp)
    elseif abs(xi) < vth/beta
        return 0.0
    else
        return Q_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Remis_murca_p(vn, vp, xi) * stepfunc(xi, vth, alpha)
    end
    
end

function Rate_murca_n(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                      xi::Float64, alpha=10.0, beta=100.0)

    vth = 3*vn + vp
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Rate_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return 0.0
    elseif abs(xi) < vth/beta
        return 0.0
    else
        return Rate_murca_n(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Rrate_murca_n(vn, vp, xi) * stepfunc(xi, vth, alpha)
    end

end

function Rate_murca_p(T::Float64, mstn::Float64, mstp::Float64, mstl::Float64, kFn::Float64, kFp::Float64, kFl::Float64,
                      SFtype_n::String, SFtype_p::String, vn::Float64, vp::Float64,
                      xi::Float64, alpha=10.0, beta=100.0)
    vth = vn + 3*vp
    if vth < threshold
        # gap size is smaller than thermal fluctuation: essentially normal fluid
        return Rate_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi)
    elseif abs(xi) < threshold
        # superfluid, but still beta-equilibrium: equilibrium one
        return 0.0
    elseif abs(xi) < vth/beta
        return 0.0
    else
        return Rate_murca_p(T, mstn, mstp, mstl, kFn, kFp, kFl, xi) * Rrate_murca_p(vn, vp, xi) * stepfunc(xi, vth, alpha)
    end
end

"""
Phase space / reduction factors
T != 0
"""

vnxis_murca_n = readdlm("../number_table/murca_n_vn_over_xi.dat", Float64, comments=true)
vpxis_murca_n = readdlm("../number_table/murca_n_vp_over_xi.dat", Float64, comments=true)
vnxis_murca_p = readdlm("../number_table/murca_p_vn_over_xi.dat", Float64, comments=true)
vpxis_murca_p = readdlm("../number_table/murca_p_vp_over_xi.dat", Float64, comments=true)

logxis = 0:0.1:2

Rrate_murca_n_nonzero_spls = [Spline2D(vnxis_murca_n[:,1], vpxis_murca_n[:,1], 
                                       readdlm("../number_table/Rrate_murca_n_SFnp_nonzeroT_logxi_$(logxi).dat", Float64, comments=true), kx=1, ky=1)
                              for logxi=logxis]

Remis_murca_n_nonzero_spls = [Spline2D(vnxis_murca_n[:,1], vpxis_murca_n[:,1], 
                                       readdlm("../number_table/Remis_murca_n_SFnp_nonzeroT_logxi_$(logxi).dat", Float64, comments=true), kx=1, ky=1)
                              for logxi=logxis]

Rrate_murca_p_nonzero_spls = [Spline2D(vnxis_murca_p[:,1], vpxis_murca_p[:,1], 
                                       readdlm("../number_table/Rrate_murca_p_SFnp_nonzeroT_logxi_$(logxi).dat", Float64, comments=true), kx=1, ky=1)
                              for logxi=logxis]

Remis_murca_p_nonzero_spls = [Spline2D(vnxis_murca_p[:,1], vpxis_murca_p[:,1], 
                                       readdlm("../number_table/Remis_murca_p_SFnp_nonzeroT_logxi_$(logxi).dat", Float64, comments=true), kx=1, ky=1)
                              for logxi=logxis]

function Rrate_murca_n_nonzero_intp(vn::Float64, vp::Float64, xi::Float64)
    logxi = log10(xi)
    idx1 = (logxis .> (logxi-0.1)) .& (logxis .<= logxi)
    idx2 = (logxis .< (logxi+0.1)) .& (logxis .>= logxi)

    vn_over_xi = vn/xi
    vp_over_xi = vp/xi

    if logxi<min(logxis...) || logxi>max(logxis...) || vn_over_xi>max(vnxis_murca_n...) || vp_over_xi>max(vpxis_murca_n...)
        return 0.0
    
    elseif idx1 == idx2
        return Rrate_murca_n_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        
    else 
        logxi1 = logxis[idx1][1]
        logxi2 = logxis[idx2][1]

        R1 = Rrate_murca_n_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        R2 = Rrate_murca_n_nonzero_spls[idx2][1](vn_over_xi, vp_over_xi)
        
        #return R1 + (R2-R1)/(logxi2-logxi1) * (logxi - logxi1)
        return R1 + (R2-R1)/(exp10(logxi2)-exp10(logxi1)) * (xi - exp10(logxi1))
    end
end

function Remis_murca_n_nonzero_intp(vn::Float64, vp::Float64, xi::Float64)
    logxi = log10(xi)
    idx1 = (logxis .> (logxi-0.1)) .& (logxis .<= logxi)
    idx2 = (logxis .< (logxi+0.1)) .& (logxis .>= logxi)

    vn_over_xi = vn/xi
    vp_over_xi = vp/xi

    if logxi<min(logxis...) || logxi>max(logxis...) || vn_over_xi>max(vnxis_murca_n...) || vp_over_xi>max(vpxis_murca_n...)
        return 0.0
    
    elseif idx1 == idx2
        return Remis_murca_n_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        
    else 
        logxi1 = logxis[idx1][1]
        logxi2 = logxis[idx2][1]

        R1 = Remis_murca_n_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        R2 = Remis_murca_n_nonzero_spls[idx2][1](vn_over_xi, vp_over_xi)
        
        #return R1 + (R2-R1)/(logxi2-logxi1) * (logxi - logxi1)
        return R1 + (R2-R1)/(exp10(logxi2)-exp10(logxi1)) * (xi - exp10(logxi1))
    end
end

function Rrate_murca_p_nonzero_intp(vn::Float64, vp::Float64, xi::Float64)
    logxi = log10(xi)
    idx1 = (logxis .> (logxi-0.1)) .& (logxis .<= logxi)
    idx2 = (logxis .< (logxi+0.1)) .& (logxis .>= logxi)

    vn_over_xi = vn/xi
    vp_over_xi = vp/xi

    if logxi<min(logxis...) || logxi>max(logxis...) || vn_over_xi>max(vnxis_murca_p...) || vp_over_xi>max(vpxis_murca_p...)
        return 0.0
    
    elseif idx1 == idx2
        return Rrate_murca_p_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        
    else 
        logxi1 = logxis[idx1][1]
        logxi2 = logxis[idx2][1]

        R1 = Rrate_murca_p_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        R2 = Rrate_murca_p_nonzero_spls[idx2][1](vn_over_xi, vp_over_xi)
        
        #return R1 + (R2-R1)/(logxi2-logxi1) * (logxi - logxi1)
        return R1 + (R2-R1)/(exp10(logxi2)-exp10(logxi1)) * (xi - exp10(logxi1))
    end
end

function Remis_murca_p_nonzero_intp(vn::Float64, vp::Float64, xi::Float64)
    logxi = log10(xi)
    idx1 = (logxis .> (logxi-0.1)) .& (logxis .<= logxi)
    idx2 = (logxis .< (logxi+0.1)) .& (logxis .>= logxi)

    vn_over_xi = vn/xi
    vp_over_xi = vp/xi

    if logxi<min(logxis...) || logxi>max(logxis...) || vn_over_xi>max(vnxis_murca_p...) || vp_over_xi>max(vpxis_murca_p...)
        return 0.0
    
    elseif idx1 == idx2
        return Remis_murca_p_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        
    else 
        logxi1 = logxis[idx1][1]
        logxi2 = logxis[idx2][1]

        R1 = Remis_murca_p_nonzero_spls[idx1][1](vn_over_xi, vp_over_xi)
        R2 = Remis_murca_p_nonzero_spls[idx2][1](vn_over_xi, vp_over_xi)
        
        #return R1 + (R2-R1)/(logxi2-logxi1) * (logxi - logxi1)
        return R1 + (R2-R1)/(exp10(logxi2)-exp10(logxi1)) * (xi - exp10(logxi1))
    end
end


"""
Phase space / reduction factors
T = 0 approx
!! xarray and yarray must be flipped for proton branch!!
"""

Remis_murca_n_table = readdlm("../number_table/Remis_murca_n.dat", Float64, comments=true)
Remis_murca_n_xarr_table = readdlm("../number_table/Remis_murca_n_xarray.dat", Float64, comments=true)
Remis_murca_n_yarr_table = readdlm("../number_table/Remis_murca_n_yarray.dat", Float64, comments=true)

Remis_murca_p_table = readdlm("../number_table/Remis_murca_p.dat", Float64, comments=true)
Remis_murca_p_xarr_table = readdlm("../number_table/Remis_murca_p_yarray.dat", Float64, comments=true)
Remis_murca_p_yarr_table = readdlm("../number_table/Remis_murca_p_xarray.dat", Float64, comments=true)

Rrate_murca_n_table = readdlm("../number_table/Rrate_murca_n.dat", Float64, comments=true)
Rrate_murca_n_xarr_table = readdlm("../number_table/Rrate_murca_n_xarray.dat", Float64, comments=true)
Rrate_murca_n_yarr_table = readdlm("../number_table/Rrate_murca_n_yarray.dat", Float64, comments=true)

Rrate_murca_p_table = readdlm("../number_table/Rrate_murca_p.dat", Float64, comments=true)
Rrate_murca_p_xarr_table = readdlm("../number_table/Rrate_murca_p_yarray.dat", Float64, comments=true)
Rrate_murca_p_yarr_table = readdlm("../number_table/Rrate_murca_p_xarray.dat", Float64, comments=true)

Remis_murca_n_spl = Spline2D(Remis_murca_n_xarr_table[1,:], Remis_murca_n_yarr_table[1,:], Remis_murca_n_table, kx=1, ky=1)
Remis_murca_p_spl = Spline2D(Remis_murca_p_xarr_table[1,:], Remis_murca_p_yarr_table[1,:], Remis_murca_p_table, kx=1, ky=1)
Rrate_murca_n_spl = Spline2D(Rrate_murca_n_xarr_table[1,:], Rrate_murca_n_yarr_table[1,:], Rrate_murca_n_table, kx=1, ky=1)
Rrate_murca_p_spl = Spline2D(Rrate_murca_p_xarr_table[1,:], Rrate_murca_p_yarr_table[1,:], Rrate_murca_p_table, kx=1, ky=1)

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

