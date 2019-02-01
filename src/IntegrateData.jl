module IntegrateData

"""
Simpson integration
"""

export integrate_data

function integrate_data(x::Array{Float64,1}, y::Array{Float64,1})
    """
    Assume
    - length(x) == length(y)
    - x is ascending
    - x[i+1] - x[i] is same for any i
    
    If length(x) is even (odd interval), we use trapezoidal rule for the first and last two, and take average
    """
    
    N = length(x)
    h = x[2] - x[1]
    if N % 2 == 1
        return h/3.0 * (y[1] + y[N] + 2*sum(y[3:2:N-2]) + 4*sum(y[2:2:N-1]))
    else
        integ_first = h/3.0 * (y[1] + y[N-1] + 2*sum(y[3:2:N-3]) + 4*sum(y[2:2:N-3])) + 0.5*h*(y[N-1]+y[N])
        integ_last = 0.5*h*(y[1]+y[2]) + h/3.0 * (y[2] + y[N] + 2*sum(y[4:2:N-2]) + 4*sum(y[3:2:N-1]))
        return (integ_first + integ_last) / 2
    end
end

function integrate_data(x::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, y::Array{Float64,1})
    """
    Assume
    - length(x) == length(y)
    - x is ascending
    - x[i+1] - x[i] is same for any i
    
    If length(x) is even (odd interval), we use trapezoidal rule for the first and last two, and take average
    """
    
    N = length(x)
    h = x[2] - x[1]
    if N % 2 == 1
        return h/3.0 * (y[1] + y[N] + 2*sum(y[3:2:N-2]) + 4*sum(y[2:2:N-1]))
    else
        integ_first = h/3.0 * (y[1] + y[N-1] + 2*sum(y[3:2:N-3]) + 4*sum(y[2:2:N-3])) + 0.5*h*(y[N-1]+y[N])
        integ_last = 0.5*h*(y[1]+y[2]) + h/3.0 * (y[2] + y[N] + 2*sum(y[4:2:N-2]) + 4*sum(y[3:2:N-1]))
        return (integ_first + integ_last) / 2
    end
end

end
