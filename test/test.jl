push!(LOAD_PATH, "../src/")

using Setup
using ODESolvers

function main()
    println(PROGRAM_FILE," start!!")

    model, core, env, var = setup("../src/sample2.ini")
    
    sol = cooling(model, core, env, var, (0.1, 2e6))
    @show sol
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

