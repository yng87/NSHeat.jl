"""
Functions to initialize the variables from ini file.
"""
# Input parameters by hand.
function setup(modelname::AbstractString, eos::AbstractString, tov::AbstractString, dMoverM::Float64, del_slice::Float64,
               Tinf0::Float64, tyr0::Float64, eta_e_inf0::Float64, eta_mu_inf0::Float64,
               SFtype_n::AbstractString, gapmodel_n::AbstractString, SFtype_p::AbstractString, gapmodel_p::AbstractString,
               noneq::Bool, P0::Float64, Pnow::Float64, Pdotnow::Float64,
               solver::AbstractString, tyrf::Float64, reltol::Float64, abstol::Float64, dt::Float64,
               output_dir::AbstractString,
               DM_heating::Bool=false,
               ann_fraction::Float64=1.0, f_capture::Float64=1.0, v_DM::Float64=230.0, rho_DM::Float64=0.42)

    model = ModelParams(modelname, eos, tov, dMoverM, del_slice,
                        SFtype_n, SFtype_p, gapmodel_n, gapmodel_p,
                        noneq, Pnow, Pdotnow, P0,
                        solver, tyrf, reltol, abstol, dt,
                        output_dir,
                        DM_heating,
                        ann_fraction, f_capture, v_DM, rho_DM)

    core = set_core_params(model)
    env = set_envelope(model)
    
    var = StarVariables(tyr0, Tinf0, eta_e_inf0, eta_mu_inf0)

    set_Tlocal(core, var)
    var.vn = similar(var.Tlocal)
    var.vp = similar(var.Tlocal)

    return model, core, env, var

end

# Input parameters by hand.
# The same as above except it receives keyword arguments
function setup(;modelname::AbstractString, eos::AbstractString, tov::AbstractString, dMoverM::Float64, del_slice::Float64,
               Tinf0::Float64, tyr0::Float64, eta_e_inf0::Float64, eta_mu_inf0::Float64,
               SFtype_n::AbstractString, gapmodel_n::AbstractString, SFtype_p::AbstractString, gapmodel_p::AbstractString,
               noneq::Bool, P0::Float64, Pnow::Float64, Pdotnow::Float64,
               solver::AbstractString, tyrf::Float64, reltol::Float64, abstol::Float64, dt::Float64,
               output_dir::AbstractString,
               DM_heating::Bool=false,
               ann_fraction::Float64=1.0, f_capture::Float64=1.0, v_DM::Float64=230.0, rho_DM::Float64=0.42)

    model = ModelParams(modelname, eos, tov, dMoverM, del_slice,
                        SFtype_n, SFtype_p, gapmodel_n, gapmodel_p,
                        noneq, Pnow, Pdotnow, P0,
                        solver, tyrf, reltol, abstol, dt,
                        output_dir,
                        DM_heating,
                        ann_fraction, f_capture, v_DM, rho_DM)

    core = set_core_params(model)
    env = set_envelope(model)
    
    var = StarVariables(tyr0, Tinf0, eta_e_inf0, eta_mu_inf0)

    set_Tlocal(core, var)
    var.vn = similar(var.Tlocal)
    var.vp = similar(var.Tlocal)

    return model, core, env, var

end

# input parameters by ini card.
function setup(filename::String)
    conf = ConfParse(filename) # ini is recommended
    parse_conf!(conf)
    
    modelname = retrieve(conf, "profile", "modelname")
    # starmodel
    eos = retrieve(conf, "starmodel", "eos")
    tov = retrieve(conf, "starmodel", "tov")
    dMoverM = retrieve(conf, "starmodel", "dMoverM", Float64)
    del_slice = retrieve(conf, "starmodel", "del_slice", Float64)
    # initial condition
    Tinf0 = retrieve(conf, "initial condition", "Tinf0", Float64)
    tyr0 = retrieve(conf, "initial condition", "tyr0", Float64)
    eta_e_inf0 = retrieve(conf, "initial condition", "eta_e_inf0", Float64)
    eta_mu_inf0 = retrieve(conf, "initial condition", "eta_mu_inf0", Float64)
    # neutron superfluidity
    SFtype_n = retrieve(conf, "neutron", "type")
    gapmodel_n = retrieve(conf, "neutron", "gap")
    # proton superfluidity
    SFtype_p = retrieve(conf, "proton", "type")
    gapmodel_p = retrieve(conf, "proton", "gap")
    # Rotochemical heating
    noneq = retrieve(conf, "rotochemical", "noneq", Bool)
    P0 = retrieve(conf, "rotochemical", "P0", Float64)
    Pnow = retrieve(conf, "rotochemical", "Pnow", Float64)
    Pdotnow = retrieve(conf, "rotochemical", "Pdotnow", Float64)
    # solver
    solver = retrieve(conf, "ODE", "solver")
    tyrf = retrieve(conf, "ODE", "tyrf", Float64)
    reltol = retrieve(conf, "ODE", "reltol", Float64)
    abstol = retrieve(conf, "ODE", "abstol", Float64)
    dt = retrieve(conf, "ODE", "dt", Float64)
    # output
    output_dir = retrieve(conf, "output", "output_dir")

    # For optional block
    if lowercase("DM") in keys(conf._data)
        DM_heating = retrieve(conf, "DM", "DM_heating", Bool)
        ann_fraction = retrieve(conf, "DM", "ann_fraction", Float64)
        f_capture = retrieve(conf, "DM", "f_capture", Float64)
        v_DM = retrieve(conf, "DM", "v_DM", Float64)
        rho_DM = retrieve(conf, "DM", "rho_DM", Float64)

    else
        DM_heating=false
        ann_fraction = 1.0
        f_capture = 1.0
        v_DM = 230.0
        rho_DM = 0.42
    end
    
    model, core, env, var = setup(modelname, eos, tov, dMoverM, del_slice,
                                  Tinf0, tyr0, eta_e_inf0, eta_mu_inf0,
                                  SFtype_n, gapmodel_n, SFtype_p, gapmodel_p,
                                  noneq, P0, Pnow, Pdotnow,
                                  solver, tyrf, reltol, abstol, dt,
                                  output_dir,
                                  DM_heating,
                                  ann_fraction, f_capture, v_DM, rho_DM)

    return model, core, env, var
end

function read_eos_core(eos_core::String)
    path_eos_core = nsheat_path * "/EOS_data/" * eos_core
    #"Rho, Press, nbar, Ye, Ymu, Yn, Yp, Yla, Ysm, Ys0, Ysp, mstp, mstn, mstla, mstsm, msts0, mstsp"
    return readdlm(path_eos_core, comments=true, comment_char='#') 
end

function read_tov(tov::String)
    path_tov = nsheat_path * "/TOV_data/Profile/" * tov
    return readdlm(path_tov, skipstart=7)
end

function get_W_Z(M::Float64, csv_file::String)
    table = readdlm(nsheat_path * "/number_table/" * csv_file, ',',  Float64, comments=true)
    sorted_idx = sortperm(table[:,1])
    x_sorted = table[:,1][sorted_idx]
    y_sorted = table[:,2][sorted_idx]
    y_spl1 = Spline1D(x_sorted, y_sorted, k=1, bc="extrapolate")
    # Note the return does not contain multiplicative factor such as 1e-61
    return y_spl1(M)
end

function set_core_params(model::ModelParams)
    eos = read_eos_core(model.EOS)
    tov = read_tov(model.TOV)
    # Get quantity from EOS, as a function of baryon density nB.
    # Reverse is needed for x to be ascending
    #Ye_intp = interpolate(eos[:,4], BSpline(Linear()))
    Ye_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,4]))
    Ymu_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,5]), k=1) # k>=2 does not work since Ymu = 0 near core-crust boundary
    Yn_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,6]))
    Yp_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,7]))
    mstp_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,12]))
    mstn_spl = Spline1D(reverse(eos[:,3]), reverse(eos[:,13]))
    

    # Get quantity from TOV table, as a function of radius
    # First, get only core params from TOV
    r_data = tov[:,2]
    nB_min_in_core = min(eos[:,3]...)
    idx_tov_core = tov[:,3] .> nB_min_in_core
    r_core_data = tov[:,2][idx_tov_core]

    # User defined slices
    del_slice = model.del_slice
    r_core = r_core_data[1]:del_slice:r_core_data[end]
        
    # Then, interpolate 
    nB_spl = Spline1D(r_core_data, tov[:,3][idx_tov_core]) #Baryon density 1/fm^3
    ephi_spl = Spline1D(r_core_data, exp.(tov[:,7][idx_tov_core])) #exp(phi)
    encl_mass_spl = Spline1D(r_core_data, tov[:,6][idx_tov_core])

    # Compute particle quantities
    nB_arr = nB_spl.(r_core)
    nn_arr = nB_arr .* Yn_spl.(nB_arr)
    np_arr = nB_arr .* Yp_spl.(nB_arr)
    ne_arr = nB_arr .* Ye_spl.(nB_arr)
    nmu_arr = nB_arr .* Ymu_spl.(nB_arr)

    kFn_arr = ( (3*pi^2) .* nn_arr ).^(1.0/3.0)
    kFp_arr = ( (3*pi^2) .* np_arr ).^(1.0/3.0)
    kFe_arr = ( (3*pi^2) .* ne_arr ).^(1.0/3.0)
    kFmu_arr = ( (3*pi^2) .* nmu_arr ).^(1.0/3.0)

    Tc_n = set_Tc_n(model, kFn_arr)
    Tc_p = set_Tc_p(model, kFp_arr)

    # Interpolate Z and W
    tov = read_tov(model.TOV)
    M = tov[:,6][end]
    Znpe = get_W_Z(M, "Znpe.csv") * 1e-61
    Znpmu = get_W_Z(M, "Znpmu.csv") * 1e-61
    Znp = get_W_Z(M, "Znp.csv") * 1e-61
    Wnpe = get_W_Z(M, "Wnpe.csv") * 1e-13
    Wnpmu = get_W_Z(M, "Wnpmu.csv") * 1e-13
    # Intialize StarCoreParams
    # Note:
    # - mst convention is different from NSCool
    #   mst is effective mass itself, not ratio
    # - radius unit is changed from [m] to [cm]
    # - we add regulator 1e-10 for volume element, which does not affect r != 0
    core = StarCoreParams(r_core .* 1e2,
                          ephi_spl.(r_core),
                          (4*pi) .* r_core.^2 .* sqrt.(r_core ./ (r_core .- (2*G).*encl_mass_spl.(r_core).*Msun./c^2 .+ 1e-10)) .* 1e4,
                          nB_arr,
                          mstn_spl.(nB_arr) .* mn,
                          mstp_spl.(nB_arr) .* mp,
                          sqrt.(me^2 .+ kFe_arr.^2 .* hbarc^2),
                          sqrt.(mmu^2 .+ kFmu_arr.^2 .* hbarc^2),
                          kFn_arr,
                          kFp_arr,
                          kFe_arr,
                          kFmu_arr,
                          nn_arr,
                          np_arr,
                          ne_arr,
                          nmu_arr,
                          Tc_n,
                          Tc_p,
                          Znpe,
                          Znpmu,
                          Znp,
                          Wnpe,
                          Wnpmu)

    return core
end

function set_envelope(model::ModelParams)
    # surface quantity
    tov = read_tov(model.TOV)
    ephi_surface = exp(tov[:,7][end])
    M = tov[:,6][end] * Msun
    R = tov[:,2][end]
    g_surface = G*M/R^2 * sqrt(R / (R-2*G*M/c^2)) * 1e2 #[cm s^-2]

    env = EnvelopeParams(ephi_surface, g_surface, R)
    return env
end

