"""
We use four structs to manage calculation efficiently.

ModelParams: parameters which describe the NS model and some other option for solving ODE and input/output.
EnvelopeParams: parameters releveant to surface-core temperature relation, and photon emission
StarCoreParams: parameters which governs the NS core structure at T=0.
                Once, fixed by constructor, these are immutable.
StarVariables: parameters for temperature evolution. 
               Only these params are mutable and change by ODE.
"""

struct ModelParams
    # model name
    modelname::String
    # Neutron star model parameter
    EOS::String
    TOV::String
    dMoverM::Float64
    del_slice::Float64
    SFtype_n::String
    SFtype_p::String
    gapmodel_n::String
    gapmodel_p::String
    noneq::Bool
    Pnow::Float64
    Pdotnow::Float64
    P0::Float64
    Znpe::Float64
    Znpmu::Float64
    Znp::Float64
    Wnpe::Float64
    Wnpmu::Float64
    # ODE
    solver::String
    tyrf::Float64
    reltol::Float64
    abstol::Float64
    dt::Float64
    # output
    output_dir::String
    # DM heating parameters
    ann_fraction::Float64 # fraction of annihilation energy converted to heating
    f_capture::Float64 # capture probability
    v_DM::Float64 # km/s
    rho_DM::Float64 # GeV/cm^3
end

struct EnvelopeParams
    ephi_surface::Float64
    g_surface::Float64
    R::Float64
end

struct StarCoreParams
    r_core::Array{Float64,1}
    ephi::Array{Float64,1}
    volume_elm::Array{Float64,1}
    nB::Array{Float64,1}
    mstn::Array{Float64,1}
    mstp::Array{Float64,1}
    mste::Array{Float64,1}
    mstmu::Array{Float64,1}
    kFn::Array{Float64,1}
    kFp::Array{Float64,1}
    kFe::Array{Float64,1}
    kFmu::Array{Float64,1}
    nn::Array{Float64,1}
    np::Array{Float64,1}
    ne::Array{Float64,1}
    nmu::Array{Float64,1}
    Tc_n::Array{Float64,1}
    Tc_p::Array{Float64,1}
end

mutable struct StarVariables
    t::Float64 #[yr]
    Tinf::Float64
    eta_e_inf::Float64 #[erg]
    eta_mu_inf::Float64 #[erg]
    # The params below are not initialized by constructor
    Tlocal::Array{Float64,1} #[K}
    vn::Array{Float64,1}
    vp::Array{Float64,1}
    Omega::Float64 #[1/s]
    Omega_dot::Float64 #[1/s^2]
    StarVariables(t, Tinf, eta_e_inf, eta_mu_inf) = new(t, Tinf, eta_e_inf, eta_mu_inf)
end

function set_Tlocal(core::StarCoreParams, var::StarVariables)
    var.Tlocal = var.Tinf ./ core.ephi
end

