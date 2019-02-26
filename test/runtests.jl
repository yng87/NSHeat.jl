push!(LOAD_PATH, "../src/")

using NSHeat
using Logging
using Test

function test()
    
    model, core, env, var = setup("sample.ini")
    sol1 = heating(model, core, env, var)
    sol2 = cooling(model, core, env, var)

    return (sol1[end]==:Success) && (sol2[end]==:Success)
end

@test test()

