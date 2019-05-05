function L_DM(model::ModelParams, core::StarCoreParams)
    # DM heating Luminosity
    # Kouvaris 2008, 10.1103/PhysRevD.77.023006
    
    NSmass = parse(Float64, split(split(model.TOV, "_")[end], ".")[1])
    vesc = 0.6425 * sqrt(NSmass/1.4) * sqrt(10/(core.r_core[end]*1e-5))
    GR_corr = sqrt(6.0/pi) # GR correction
    gamma = 1.0/sqrt(1-vesc^2)
    return 3.41e22 * (model.ann_fraction+gamma-1.0) * model.f_capture * (core.r_core[end]*1e-5/10) * (230.0/model.v_DM) * (model.rho_DM/0.42) * (NSmass/1.4) * GR_corr
end

