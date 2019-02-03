module SpinDown

export set_Omega, set_Omega_dot

push!(LOAD_PATH, "./")

using NeutronStar

function set_Omega(model::ModelParams, var::StarVariables)
    # 1/s
    var.Omega = 2*pi / sqrt(model.P0^2 + 2*model.Pnow*model.Pdotnow*var.t)
    return nothing
end

function set_Omega_dot(model::ModelParams, var::StarVariables)
    # 1/s^2
    var.Omega_dot = -2*pi*model.Pnow*model.Pdotnow / (sqrt(model.P0^2 + 2*model.Pnow*model.Pdotnow*var.t))^3
    return nothing
end

end
