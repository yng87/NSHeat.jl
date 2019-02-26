"""
Spin-down with pulsar breaking index n=3.
"""

function set_Omega(model::ModelParams, var::StarVariables)
    # 1/s
    var.Omega = 2*pi / sqrt(model.P0^2 + 2*model.Pnow*model.Pdotnow*var.t*yrTosec)
    return nothing
end

function set_Omega_dot(model::ModelParams, var::StarVariables)
    # 1/s^2
    var.Omega_dot = -2*pi*model.Pnow*model.Pdotnow / (sqrt(model.P0^2 + 2*model.Pnow*model.Pdotnow*var.t*yrTosec))^3
    return nothing
end

