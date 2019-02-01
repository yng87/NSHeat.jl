module SpecHeat

export spec_heat

push!(LOAD_PATH, "./")
include("./PhysicalConstants.jl")

# erg*MeV/fm -> fm^-3 -> cm^-3
converter = 1.0/MeVToerg * MeVTofminv^2 * 1e39

function spec_heat(mst::Float64, kF::Float64, T::Float64)
    # output is [erg/K/cm^3]
    # without converter, unit is [MeV erg^2 fm^-1 K^-1]
    return mst*kF*kB^2*T / 3 * converter
end

function spec_heat(mst::Float64, kF::Float64, T::Float64, v::Float64, SFtype::String)
    # For superfluid species
    if v == 0.0
        return spec_heat(mst, kF, T)
    elseif SFtype == "1S0"
        return spec_heat(mst, kF, T) * RA(v)
    elseif SFtype == "3P2m0"
        return spec_heat(mst, kF, T) * RB(v)
    elseif SFtype == "3P2m2"
        return spec_heat(mst, kF, T) * RC(v)
    else
        println("spec_heat: wrong superfluid type")
        return 0.0
    end
end

RA(v::Float64) = (0.4186 + sqrt(1.007^2+0.5010^2*v^2))^(2.5) * exp(1.456 - sqrt(1.456^2 + v^2))
RB(v::Float64) = (0.6893 + sqrt(0.790^2+0.2824^2*v^2))^2 * exp(1.934 - sqrt(1.934^2+v^2))
RC(v::Float64) = (2.188 - (9.537e-5)^2*v^2 + 0.1491^4*v^4)/(1.0 + 0.2846^2*v^2 + 0.01335^4*v^4 + 0.1815^6*v^6)


function main()
    println(PROGRAM_FILE," start!!")

    @show spec_heat(900., 0.1, 1e10, 0.0, "3P2m0")
    @show spec_heat(900., 0.1, 1e10, 0.00001, "3P2m0")
    println(PROGRAM_FILE," finish!!")
end

if occursin(PROGRAM_FILE, @__FILE__)
    @time main()
end

    

end
