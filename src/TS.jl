export AdaptTS

abstract type AbstractTS end

Base.@kwdef mutable struct AdaptTS <: AbstractTS
    Δt::Real = 1e-3
    errₘₐₓ::Real = 1e-3
    rel::Bool = false
end
