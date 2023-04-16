export Lorenz

Base.@kwdef struct Lorenz <: AbstractProblem
    a::Real = 16
    r::Real = 45
    b::Real = 4
end

dim(::Lorenz) = 3

function ∂ₜ(l::Lorenz, ::NoForcing)
    (; a, b, r) = l
    return function (state, _)
        x, y, z = state[1:3]
        return Float64[
            a*(y-x),
            r*x-y-x*z,
            x*y-b*z
        ]
    end
end

function Dₓ(l::Lorenz)
    (; a, b, r) = l
    return function (state, _)
        x, y, z = state[1:3]
        Float64[
            -a a 0
            r-z -1 -x
            y x -b
        ]
    end
end