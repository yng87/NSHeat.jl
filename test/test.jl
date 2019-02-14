push!(LOAD_PATH, "../src/")

using NSHeat
using Logging

function main()
    println(PROGRAM_FILE," start!!")
    
    model, core, env, var = setup(ARGS[1])
    if model.noneq == true
        println("Heating:")
        sol = heating(model, core, env, var)
        write_ini(sol, model)
        output_T(sol, model, core, env, var)
        output_LC(sol, model, core, env, var)
    else
        println("Cooling:")
        sol = cooling(model, core, env, var)
        write_ini(sol, model)
        output_T(sol, model, core, env, var)
        output_LC(sol, model, core, env, var)
    end
    
    println(PROGRAM_FILE," finish!!")

    return sol[end]
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

