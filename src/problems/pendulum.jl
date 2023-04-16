export Pendulum, DampedPendulum

Base.@kwdef struct Pendulum <: AbstractProblem
    "Mass (kg)"
    m::Real = 0.1
    "Length (m)"
    ℓ::Real = 0.1
    "Gravitational Acceleration (m/s²)"
    g::Real = 9.8
end

Base.@kwdef struct DampedPendulum <: AbstractProblem
    "Mass (kg)"
    m::Real = 0.1
    "Length (m)"
    ℓ::Real = 0.1
    "Coefficient of Friction"
    β::Real = 0
    "Gravitational Acceleration (m/s²)"
    g::Real = 9.8
end

convert(::Type{DampedPendulum}, p::Pendulum) = DampedPendulum(p.m, p.ℓ, 0, p.g)
dim(::Pendulum) = 2

function ∂ₜ(dampedPendulum::DampedPendulum, f::AbstractForcing=NoForcing())
    (; m, ℓ, β, g) = dampedPendulum
    return (x, t) -> [x[2], (force(f, t) / (m * ℓ) - β / m * x[2] - g / ℓ * sin(x[1]))]
end