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

function set_Omega(model::ModelParams, var::StarVariables, n::Int)
    var.Omega = ifelse(n==1,
                       2*pi/model.P0 * exp(-model.Pdotnow/model.Pnow*var.t*yrTosec),
                       2*pi/(model.P0^(n-1) + model.Pdotnow*model.Pnow^(n-2)*(n-1)*var.t*yrTosec)^(1/(n-1)))
    return nothing
end

function set_Omega_dot(model::ModelParams, var::StarVariables, n::Int)
    var.Omega_dot = ifelse(n==1,
                           -model.Pdotnow/model.Pnow*2*pi/model.P0 * exp(-model.Pdotnow/model.Pnow*var.t*yrTosec),
                           -2*pi*model.Pdotnow*model.Pnow^(n-2)/(model.P0^(n-1) + model.Pdotnow*model.Pnow^(n-2)*(n-1)*var.t*yrTosec)^(n/(n-1.0)))
    return nothing
end
