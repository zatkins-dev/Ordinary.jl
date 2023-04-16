export Rössler, Rossler

Base.@kwdef struct Rössler <: AbstractProblem
    a::Real
    b::Real
    c::Real
end

const Rossler = Rössler

dim(::Rössler) = 3

function ∂ₜ(r::Rössler, ::NoForcing)
    (; a, b, c) = r
    (x, _) -> [
        -(x[2] + x[3]),
        x[1] + a * x[2],
        b + x[3] * (x[1] - c)
    ]
end

function Dₓ(r::Rössler)
    (; a, b, c) = r
    return (x, _) ->
        Float64[
            0 -1 -1
            1 a 0
            x[3] 0 x[1]-c
        ]
end