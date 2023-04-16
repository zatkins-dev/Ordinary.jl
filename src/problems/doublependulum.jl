export DoublePendulum, ∂ₜ

Base.@kwdef struct DoublePendulum <: AbstractProblem
    "Mass of First Pendulum (kg)"
    m₁::Real = 0.1
    "Length of First Pendulum (m)"
    ℓ₁::Real = 0.1
    "Mass of Second Pendulum (kg)"
    m₂::Real = 0.1
    "Length of Second Pendulum (m)"
    ℓ₂::Real = 0.1
    "Gravitational Acceleration (m/s²)"
    g::Real = 9.8
end

# function C₁(dp::DoublePendulum, x)
#     θ₁, θ₂, p₁, p₂ = x
#     (p₁ * p₂ * sin(θ₁ - θ₂)) / (dp.ℓ₁ * dp.ℓ₂ * (dp.m₁ + dp.m₂ * sin(θ₁ - θ₂)^2))
# end

# function C₂(dp::DoublePendulum, x)
#     θ₁, θ₂, p₁, p₂ = x
#     (dp.ℓ₂^2 * dp.m₁^2 * p₁ + dp.ℓ₁^2 * (dp.m₁ + dp.m₂) * p₂^2 - dp.ℓ₁ * dp.ℓ₂ * dp.m₂ * p₁ * p₂ * cos(θ₁ - θ₂)) / (2 * dp.ℓ₁^2 * dp.ℓ₂^2 * (dp.m₁ + dp.m₂ * sin(θ₁ - θ₂)^2)^2) * sin(2(θ₁ - θ₂))
# end

# function ∂ₜ(dp::DoublePendulum, ::NoForcing)
#     return function (x, _)
#         θ₁, θ₂, p₁, p₂ = x
#         c₁ = C₁(dp, x)
#         c₂ = C₂(dp, x)
#         return [
#             (dp.ℓ₂ * p₁ - dp.ℓ₁ * p₂ * cos(θ₁ - θ₂)) / (dp.ℓ₁^2 * dp.ℓ₂ * (dp.m₁ + dp.m₂ * sin(θ₁ - θ₂))),
#             (dp.ℓ₁ * (dp.m₁ + dp.m₂) * p₂ - dp.ℓ₂ * dp.m₂ * p₁ * cos(θ₁ - θ₂)) / (dp.ℓ₁ * dp.ℓ₂^2 * dp.m₂ * (dp.m₁ + dp.m₂ * sin(θ₁ - θ₂))),
#             -(dp.m₁ + dp.m₂) * dp.g * dp.ℓ₁ * sin(θ₁) - c₁ + c₂,
#             -dp.m₂ * dp.g * dp.ℓ₂ * sin(θ₂) + c₁ - c₂
#         ]
#     end
# end

function ∂ₜ(dp::DoublePendulum, ::NoForcing)
    (; m₁, m₂, ℓ₁, ℓ₂, g) = dp
    return function (x, _)
        θ₁, θ₂, ω₁, ω₂ = x
        return [
            ω₁,
            ω₂,
            (-g * (2m₁ + m₂) * sin(θ₁) - m₂ * g * sin(θ₁ - 2θ₂) - 2sin(θ₁ - θ₂) * m₂ * (ω₂^2 * ℓ₂ + ω₁^2 * ℓ₁ * cos(θ₁ - θ₂))) / (ℓ₁ * (2m₁ + m₂ - m₂ * cos(2(θ₁ - θ₂)))),
            2sin(θ₁ - θ₂) * (ω₁^2 * ℓ₁ * (m₁ + m₂) + g * (m₁ + m₂) * cos(θ₁) + ω₂^2 * ℓ₂ * m₂ * cos(θ₁ - θ₂)) / (ℓ₂ * (2m₁ + m₂ - m₂ * cos(2(θ₁ - θ₂))))
        ]
    end
end