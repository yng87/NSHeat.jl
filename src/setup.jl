module Setup

export set_core_params

include("./NeutronStar.jl")
include("./PhysicalConstants.jl")

using .NeutronStar
using DelimitedFiles
using Dierckx

function read_eos_core(path_eos_core::String)
    #"Rho, Press, nbar, Ye, Ymu, Yn, Yp, Yla, Ysm, Ys0, Ysp, mstp, mstn, mstla, mstsm, msts0, mstsp"
    return readdlm(path_eos_core, comments=true, comment_char='#') 
end

function read_tov(path_tov::String)
    return readdlm(path_tov, skipstart=7)
end

function set_core_params(eos::Array{Float64,2}, tov::Array{Float64,2}, del_slice::Float64)
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
    
pp    # Intialize StarCoreParams
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
                          nmu_arr)

    return core
end

function main()
    println(PROGRAM_FILE," start!!")

    eos = read_eos_core("../EOS_data/APR_EOS_Cat_core.dat")
    tov = read_tov("../TOV_data/Profile/Prof_APR_Cat_1.4.dat")
    del_slice=100.0
    println(set_core_params(eos, tov, del_slice))
    
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

end
