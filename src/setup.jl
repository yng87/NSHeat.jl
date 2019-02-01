module Setup
"""
Functions to initialize the variables from ini file.
"""

export set_core_params, setup

include("./PhysicalConstants.jl")
push!(LOAD_PATH, "./")

using ConfParser
using NeutronStar
using SuperfluidGaps
using DelimitedFiles
using Dierckx

function setup(filename::String)
    conf = ConfParse(filename) # ini is recommended
    parse_conf!(conf)
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
    rotochemical = retrieve(conf, "rotochemical", "rc", Bool)
    P0 = retrieve(conf, "rotochemical", "P0", Float64)
    Pnow = retrieve(conf, "rotochemical", "Pnow", Float64)
    Pdotnow = retrieve(conf, "rotochemical", "Pdotnow", Float64)
    Znpe = retrieve(conf, "rotochemical", "Znpe", Float64)
    Znpmu = retrieve(conf, "rotochemical", "Znpmu", Float64)
    Znp = retrieve(conf, "rotochemical", "Znp", Float64)
    Wnpe = retrieve(conf, "rotochemical", "Wnpe", Float64)
    Wnpmu = retrieve(conf, "rotochemical", "Wnpmu", Float64)
    # output
    output_dir = retrieve(conf, "output", "output_dir")

    model = ModelParams(eos, tov, dMoverM, del_slice,
                        SFtype_n, SFtype_p, gapmodel_n, gapmodel_p,
                        rotochemical, Pnow, Pdotnow, P0,
                        Znpe, Znpmu, Znp, Wnpe, Wnpmu,
                        output_dir)

    core = set_core_params(model)
    
    var = StarVariables(tyr0, Tinf0, eta_e_inf0, eta_mu_inf0)
    var.Tlocal = Tinf0 ./ core.ephi
    var.vn = similar(var.Tlocal)
    var.vp = similar(var.Tlocal)

    return model, core, var
end

function read_eos_core(path_eos_core::String)
    #"Rho, Press, nbar, Ye, Ymu, Yn, Yp, Yla, Ysm, Ys0, Ysp, mstp, mstn, mstla, mstsm, msts0, mstsp"
    return readdlm(path_eos_core, comments=true, comment_char='#') 
end

function read_tov(path_tov::String)
    return readdlm(path_tov, skipstart=7)
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
    # Intialize StarCoreParams
    core = StarCoreParams(r_core,
                          ephi_spl.(r_core),
                          nB_arr,
                          mstn_spl.(nB_arr),
                          mstp_spl.(nB_arr),
                          sqrt.(1 .+ kFe_arr.^2 ./me^2),
                          sqrt.(1 .+ kFmu_arr.^2 ./mmu^2),
                          kFn_arr,
                          kFp_arr,
                          kFe_arr,
                          kFmu_arr,
                          nn_arr,
                          np_arr,
                          ne_arr,
                          nmu_arr,
                          Tc_n,
                          Tc_p)
    # Tc_{n,p} is set to 0

    return core
end

function main()
    println(PROGRAM_FILE," start!!")

    x, y, z = setup("./sample.ini")

    @show y.Tc_n
    @show y.Tc_p
    set_vp(x, y, z)
    set_vn(x, y, z)
    @show z.vp
    @show z.vn
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

end
