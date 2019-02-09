push!(LOAD_PATH, "../src/")

using Setup
using ODESolvers
using Output
using Logging

function run(model, core, env, var)
    @show model.modelname
        
    try
        if model.noneq == true
            sol = heating(model, core, env, var, 1e-10, 1e-10)
        else
            sol = cooling(model, core, env, var, 1e-10, 1e-10)
        end
        @show sol.retcode
        if sol.retcode == :Success
            write_ini(sol, model)
            output_T(sol, model, core, env, var)
            output_LC(sol, model, core, env, var)
        else
            @error "ODE solver failed" model.modelname
            @error "status" sol.retcode
        end
    catch err
        @show "Failed"
        @error "ODE solver failed" model.modelname
        @error "Error" err
    end
end
        
function main()
    println(PROGRAM_FILE," start!!")
    
    model, core, env, var = setup(ARGS[1])
    if model.noneq == true
        sol = heating(model, core, env, var, 1e-4, 1e-4)
    else
        sol = cooling(model, core, env, var, 1e-10, 1e-10)
    end
    #run(model, core, env, var)
    
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

