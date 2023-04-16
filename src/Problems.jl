export ∂ₜ, Dₓ

using LinearAlgebra

abstract type AbstractProblem end
abstract type AbstractForcing end

# Forcing terms
include("problems/forcing.jl")

# Problems
include("problems/pendulum.jl")
include("problems/lorenz.jl")
include("problems/rossler.jl")
include("problems/doublependulum.jl")


∂ₜ(p::AbstractProblem, f::Function) = ∂ₜ(p, FunctionDriver(f))
∂ₜ(p::AbstractProblem) = ∂ₜ(p, NoForcing())
Dₓ(p::AbstractProblem) = (_, _) -> I