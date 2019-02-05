module Output

export output_T, output_LC

push!(LOAD_PATH, "./")

using DifferentialEquations
using NeutronStar
using HeatCapacity
using NeutrinoLum
using PhotonLum
using SuperfluidGaps
using SpinDown

function output_T(sol, model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    filepath = model.output_dir * "temperature.dat"

    open(filepath, "w") do io
        println(io, "# t[yr] Teff_inf[K] Tinf[K] eta_e_inf[erg] eta_mu_inf[erg]")
        for (t, u) in zip(sol.t, sol.u)
            var.Tinf = u[1]
            Teff_inf = get_Teff_inf(model, env, var)
            println(io, "$t $(Teff_inf) $(u[1]) $(u[2]) $(u[3])")
        end
    end
end

function output_LC(sol, model::ModelParams, core::StarCoreParams, env::EnvelopeParams, var::StarVariables)
    filepath_L = model.output_dir * "luminosity.dat"
    filepath_C = model.output_dir * "capacity.dat"

    ioL = open(filepath_L, "w")
    ioC = open(filepath_C, "w")

    println(ioL, "# t[yr] Urca PBF Photon Heat")
    println(ioL, "# Luminosity unit = [erg/s]")
    println(ioC, "# t[yr] n p e mu")
    println(ioC, "# Heat Capacity unit = [erg/K]")
    for (t, u) in zip(sol.t, sol.u)
        var.t = t #yr
        var.Tinf = u[1]
        var.eta_e_inf = u[2] #erg
        var.eta_mu_inf = u[3] #erg
        set_Tlocal(core, var)
        set_vn(model, core, var)
        set_vp(model, core, var)
        set_Omega(model, var)
        set_Omega_dot(model, var)
        
        Lurca = L_murca_n_e(model, core, var, model.noneq) + L_murca_n_mu(model, core, var, model.noneq) + L_murca_p_e(model, core, var, model.noneq) + L_murca_p_mu(model, core, var, model.noneq)
        
        LPBF = 0.0
        if lowercase(model.SFtype_n) != "normal"
            LPBF += L_PBF_n(model, core, var)
        end
        if lowercase(model.SFtype_p) != "normal"
            LPBF += L_PBF_p(model, core, var)
        end
        
        Lphoton = L_photon(model, env, var)
        
        LHeat = 0.0
        if model.noneq == true
            Rate_e = Rate_volume_murca_n_e(model, core, var, model.noneq) + Rate_volume_murca_p_e(model, core, var, model.noneq)
            Rate_mu = Rate_volume_murca_n_mu(model, core, var, model.noneq) + Rate_volume_murca_p_mu(model, core, var, model.noneq)
            LHeat = var.eta_e_inf*Rate_e + var.eta_mu_inf*Rate_mu
        end

        Ce = get_Ce(model, core, var)
        Cmu = get_Cmu(model, core, var)
        Cn = get_Cn(model, core, var)
        Cp = get_Cp(model, core, var)
        
        println(ioL, "$t $(Lurca) $(LPBF) $(Lphoton) $(LHeat)")
        println(ioC, "$t $(Cn) $(Cp) $(Ce) $(Cmu)")
    end

    close(ioL)
    close(ioC)

end


            
end
