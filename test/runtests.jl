#push!(LOAD_PATH, "../src/")

using NSHeat
using Logging
using Test

function test()

    success = true

    println("Test normal fluid...")
    model, core, env, var = setup("sample.ini")
    sol = cooling(model, core, env, var)
    success = success && (sol[end]==:Success)

    model, core, env, var = setup("sample.ini")
    sol = heating(model, core, env, var)
    success = success && (sol[end]==:Success)
    @show success

    println("Test superfluid with radau...")
    model, core, env, var = setup("sample2.ini")
    sol = cooling(model, core, env, var)
    success = success && (sol[end]==:Success)

    model, core, env, var = setup("sample2.ini")
    sol = heating(model, core, env, var)
    success = success && (sol[end]==:Success)
    @show success

    println("Test superfluid with CVODE_BDF...")
    model, core, env, var = setup("sample3.ini")
    sol = cooling(model, core, env, var)
    success = success && (sol[end]==:Success)

    model, core, env, var = setup("sample3.ini")
    sol = heating(model, core, env, var)
    success = success && (sol[end]==:Success)
    @show success

    return success
end

@test test()

