using LinearAlgebra
export rk4, euler, solve

Optional{T} = Union{T,Missing}

include("TS.jl")

function euler(∂ₜx::Function, x::AbstractVector, t::Real, Δt::Real)
    x + Δt * ∂ₜx(x, t), t + Δt
end

function euler(∂ₜx::Function, DₓF::Function, x::AbstractVector, δ::AbstractMatrix, t::Real, Δt::Real)
    x + Δt * ∂ₜx(x, t), δ + Δt * DₓF(x, t) * δ, t + Δt
end

function rk4(∂ₜx::Function, x::AbstractVector, t::Real, Δt::Real)
    k₁ = ∂ₜx(x, t)
    k₂ = ∂ₜx(x + Δt / 2 * k₁, t + Δt / 2)
    k₃ = ∂ₜx(x + Δt / 2 * k₂, t + Δt / 2)
    k₄ = ∂ₜx(x + Δt * k₃, t + Δt)
    x + Δt * (k₁ + 2 * k₂ + 2 * k₃ + k₄) / 6, t + Δt
end

function rk4(∂ₜx::Function, DₓF::Function, x::AbstractVector, δ::AbstractMatrix{<:AbstractFloat}, t::Real, Δt::Real)
    k₁, d₁ = ∂ₜx(x, t), DₓF(x, t) * δ
    k₂, d₂ = ∂ₜx(x + 0.5Δt * k₁, t + 0.5Δt), DₓF(x + 0.5Δt * k₁, t + 0.5Δt) * (δ + 0.5Δt * d₁)
    k₃, d₃ = ∂ₜx(x + 0.5Δt * k₂, t + 0.5Δt), DₓF(x + 0.5Δt * k₂, t + 0.5Δt) * (δ + 0.5Δt * d₂)
    k₄, d₄ = ∂ₜx(x + Δt * k₃, t + Δt), DₓF(x + Δt * k₃, t + Δt) * (δ + Δt * d₃)
    x + Δt * (k₁ + 2 * k₂ + 2 * k₃ + k₄) / 6, δ + Δt * (d₁ + 2 * d₂ + 2 * d₃ + d₄) / 6, t + Δt
end

function rk4(∂ₜx::Function, x::AbstractVector, t::Real, ts::AdaptTS)
    @debug "Computing rk4 step with adaptive timestepping"
    errf(u, v) = norm(u - v, Inf) / (ts.rel ? norm(u) : 1)
    # Verification
    xₙ₊₁, tₙ₊₁ = rk4(∂ₜx, x, t, ts.Δt)
    x⁺, _ = rk4(∂ₜx, x, t, ts.Δt / 2)
    x₂, _ = rk4(∂ₜx, x⁺, t + ts.Δt / 2, ts.Δt / 2)
    err = errf(x₂, xₙ₊₁)
    @debug err
    if err ≤ ts.errₘₐₓ # Can relax Δt
        x₂, _ = rk4(∂ₜx, xₙ₊₁, t + ts.Δt, ts.Δt)
        x₁, t₁ = rk4(∂ₜx, x, t, 2ts.Δt)
        err = errf(x₂, x₁)
        while err ≤ ts.errₘₐₓ
            xₙ₊₁, tₙ₊₁ = x₁, t₁
            ts.Δt = 2ts.Δt
            @debug "[t=$(t)] Relaxing Δt" ts.Δt
            x₂, _ = rk4(∂ₜx, xₙ₊₁, t + ts.Δt, ts.Δt)
            x₁, t₁ = rk4(∂ₜx, x, t, 2ts.Δt)
            err = errf(x₂, x₁)
        end
    else # Need smaller Δt
        while err > ts.errₘₐₓ
            ts.Δt = ts.Δt / 2
            @debug "[t=$(t)] Shrinking Δt" ts.Δt
            xₙ₊₁, tₙ₊₁ = rk4(∂ₜx, x, t, ts.Δt)
            # Verification - 2 Half Steps
            x⁺, _ = rk4(∂ₜx, x, t, ts.Δt / 2)
            x₂, _ = rk4(∂ₜx, x⁺, t + ts.Δt / 2, ts.Δt / 2)
            err = errf(x₂, xₙ₊₁)
        end
    end
    xₙ₊₁, tₙ₊₁
end

function solve(solver::Function, ∂ₜx::Function; x₀::AbstractVector, t₀=0.0, Δt::Union{Real,AbstractTS}=1e-3, steps::Optional{Int}=missing, tₑ::Optional{Real}=missing, returnvectors=false)
    if !(ismissing(steps) ⊻ ismissing(tₑ))
        error("Must provide exactly one of number of steps (steps) or ending time (tₑ) as kwarg")
    end
    X = Vector{Float64}[]
    push!(X, x₀)
    V = Vector{Float64}[]
    t = Real[]
    push!(t, t₀)
    check = ismissing(tₑ) ? (ts) -> length(ts) <= steps : (ts) -> ts[end] < tₑ
    while check(t)
        xₙ₊₁, tₙ₊₁ = solver(∂ₜx, X[end], t[end], Δt)
        push!(t, tₙ₊₁)
        push!(X, xₙ₊₁)
        if returnvectors
            push!(V, X[end] - X[end-1])
        end
    end
    Xr = reduce(hcat, X)
    if returnvectors
        push!(V, zero(V[end]))
        Vr = reduce(hcat, V)
        return Xr, Vr, t
    else
        return Xr, t
    end
end

function solve(solver::Function, ∂ₜx::Function, DₓF::Function; x₀::AbstractVector, δ₀::AbstractArray=Diagonal{Float64}(1.0I, 3), t₀=0.0, Δt::Union{Real,AbstractTS}=1e-3, steps::Optional{Int}=missing, tₑ::Optional{Real}=missing, returnvectors=false)
    if !(ismissing(steps) ⊻ ismissing(tₑ))
        error("Must provide exactly one of number of steps (steps) or ending time (tₑ) as kwarg")
    end
    n = length(x₀)
    if length(δ₀) != n^2
        error("Dimension of initial variations should be n²")
    end
    if ndims(δ₀) != 2
        δ₀ = reshape(δ₀, (n, n))
    end
    X = Vector{Float64}[]
    δ = Matrix{Float64}[]
    push!(X, x₀)
    push!(δ, δ₀)
    V = Vector{Float64}[]
    t = Real[]
    push!(t, t₀)
    check = ismissing(tₑ) ? (ts) -> length(ts) <= steps : (ts) -> ts[end] < tₑ
    while check(t)
        xₙ₊₁, δₙ₊₁, tₙ₊₁ = solver(∂ₜx, DₓF, X[end], δ[end], t[end], Δt)
        push!(t, tₙ₊₁)
        push!(X, xₙ₊₁)
        push!(δ, δₙ₊₁)
        if returnvectors
            push!(V, X[end] - X[end-1])
        end
    end
    Xr = reduce(hcat, X)
    if returnvectors
        push!(V, zero(V[end]))
        Vr = reduce(hcat, V)
        return Xr, δ, Vr, t
    else
        return Xr, δ, t
    end
end