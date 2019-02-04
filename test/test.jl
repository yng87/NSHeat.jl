push!(LOAD_PATH, "../src/")

using Setup
using ODESolvers
using Plots

function main()
    println(PROGRAM_FILE," start!!")

    model, core, env, var = setup("../src/sample.ini")
    
    #sol = cooling(model, core, env, var, (0.1, 2e6))
    sol = heating(model, core, env, var, (1., 1e8), 1e-14, 1e-14)
    @show sol
    #T = map(u->u[1], sol.u)

    #pyplot()
    #plot(sol.t, T, scale=:log10)
    #show()
    #savefig("~/Desktop/test.pdf")
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

