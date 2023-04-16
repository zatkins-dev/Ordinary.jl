export CosineDriver, FunctionDriver, NoForcing, force


# Forcing structs
Base.@kwdef struct CosineDriver <: AbstractForcing
    α::Real
    A::Real
end

struct FunctionDriver <: AbstractForcing
    f::Function
end

struct NoForcing <: AbstractForcing end

# Get forcing functions
force(c::CosineDriver, t::Real) = c.A * cos(c.α * t)
force(f::Function, t::Real) = f(t)
force(::NoForcing, ::Real) = 0