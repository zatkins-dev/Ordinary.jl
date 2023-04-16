export Wolf, NearestNeighbor, DirectionalNearestNeighbor
using LinearAlgebra

abstract type AbstractNearestNeighbor end

Base.@kwdef struct Wolf
    ϵ::Float64 = 0.1
    theiler::Float64 = 0
end

Base.@kwdef mutable struct WolfResults
    w::Wolf
    NΔt::Float64 = 1
    Ls::Vector{Float64} = []
    Le::Vector{Float64} = []
    λ₁::Float64 = 0
end

Base.@kwdef mutable struct DirectionalNearestNeighbor <: AbstractNearestNeighbor
    θ::Float64 = π / 9
    max_angle::Float64 = π / 2
    weight_distance = 0.5
    bidirectional = false
end

struct NearestNeighbor <: AbstractNearestNeighbor end

function (::NearestNeighbor)(x, traj, mask, ϵ; verbose=0, kwargs...)
    dist = [norm(xk - x) for xk in traj]
    dist_mask = dist .< ϵ .&& mask
    if !any(dist_mask)
        error("No points within ϵ = $(ϵ) of point $(x'), run with higher ϵ")
    end
    min_k_in_sphere = argmin(dist[dist_mask])
    k = (1:length(dist_mask))[dist_mask][min_k_in_sphere]
    return k
end

function (dnn::DirectionalNearestNeighbor)(x, traj, mask, ϵ; last_direction::Union{Nothing,Vector}=nothing, verbose=0)
    k = -1
    if last_direction === nothing || any(isnan.(last_direction))
        k = NearestNeighbor()(x, traj, mask, ϵ; verbose=verbose)
    else
        normalize!(last_direction)
        dist = [norm(xk - x) for xk in traj]
        dist_mask = dist .< ϵ .&& mask
        if !any(dist_mask)
            error("No points within ϵ = $(ϵ) of point $(x'), run with higher ϵ")
        end
        f = dnn.bidirectional ? abs : identity
        angle = [acos(min(f((xk - x) / d ⋅ last_direction), 1)) for (xk, d) in zip(traj, dist)]
        θ = dnn.θ
        cone = BitArray(dist_mask .&& angle .< θ)
        while !any(cone) && θ < dnn.max_angle
            θ += 0.01
            @. cone = dist_mask && angle .< θ
        end
        if !any(cone)
            error("No points within ϵ = $(ϵ), θ = $(θ/π)π cone from point $(x'), direction $(last_direction'), run with higher ϵ")
        end
        if count(cone) < 10
            k = (1:length(cone))[cone][argmin(dist[cone])]
        else
            ω = dnn.weight_distance
            max_k_in_cone = argmin(ω * dist[cone] + (1 - ω) * angle[cone])
            k = (1:length(cone))[cone][max_k_in_cone]
            if verbose > 1
                @show k, dist[cone][max_k_in_cone], minimum(dist[cone]), angle[cone][max_k_in_cone], minimum(angle[cone])
                @show extrema(dist[cone]), extrema(angle[cone])
            end
        end
    end
    if verbose > 1
        @show k
    end
    return k
end

function (w::Wolf)(x::Vector, t::Vector, nn::AbstractNearestNeighbor; verbose=0)
    start = round(Int, 3 * length(t) / 8)
    results = WolfResults(w=w)
    indices = 1:length(t)
    n = start
    k, j = 0, 0
    while n < length(t)
        mask = abs.(indices .- n) .> w.theiler
        mask[end] = false
        k = nn(x[n], x, mask, w.ϵ; last_direction=k == 0 ? nothing : normalize(x[k] - x[n]), verbose=verbose)
        Ls = norm(x[n] - x[k])
        Le = Ls
        j = 0
        while Le < w.ϵ && k < length(t) && n < length(t)
            n += 1
            k += 1
            j += 1
            Le = norm(x[k] - x[n])
        end
        if Le < w.ϵ
            if n == length(t)
                n -= j
                break
            end
            continue
        end
        results.λ₁ += abs(Le - Ls) > eps() ? log2(Le) - log2(Ls) : 0
        push!(results.Ls, Ls)
        push!(results.Le, Le)
        if verbose > 0
            @show n, k, Ls, Le
        end
    end
    results.λ₁ /= max(t[n] - t[start], 1)
    return results
end

