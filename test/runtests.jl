#push!(LOAD_PATH, "../src/")

using NSHeat
using Logging
using Test

function test()
    
    model, core, env, var = setup("sample.ini")
    sol1 = heating(model, core, env, var)
    @show sol1[1][end], sol1[2][end], sol1[end]
    output_LC(sol1, model, core, env, var)
    model, core, env, var = setup("sample.ini")
    sol2 = cooling(model, core, env, var)
    @show sol2[1][end], sol2[2][end], sol2[end]

    return (sol1[end]==:Success) && (sol2[end]==:Success)
end

@test test()

