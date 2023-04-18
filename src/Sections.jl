using LinearAlgebra

export TemporalSection, Hyperplane, SpacialSection, section

abstract type AbstractSection end

struct TemporalSection <: AbstractSection
    "Period of the section"
    T::Real
end

Base.@kwdef struct Hyperplane
    "Point on plane"
    x::Vector{Real}
    "Normal vector to plane"
    n::Vector{Real}
end

function insideoutside(plane::Hyperplane, x)
    (x - plane.x) ⋅ plane.n
end

Hyperplane(x₀::Vector{Real}, F::AbstractProblem; t=0) = Hyperplane(x=x₀, n=∂ₜ(F)(x₀, t))

struct SpacialSection <: AbstractSection
    "Plane of spacial section"
    plane::Hyperplane
end

function insideoutside(Σ::TemporalSection, tₚ, t)
    (floor(Int, t / Σ.T) - floor(Int, tₚ / Σ.T)) > 0
end

function insideoutside(Σ::SpacialSection, xₚ, x)
    insideoutside(Σ.plane, xₚ) < 0 && insideoutside(Σ.plane, x) > 0
end

function section(Σ::TemporalSection, xs::Matrix, ts::Vector; interpolate::Bool=false, returnindex::Bool=false)
    indices = findall(t -> insideoutside(Σ, t[1], t[2]), collect(zip(ts[1:end-1], ts[2:end]))) .+ 1
    X = xs[:, indices]
    t = ts[indices]
    if interpolate
        dt = ts[indices] - ts[indices.-1]
        # temporaldist = (Σ.T .- mod1.(ts[indices.-1], Σ.T))
        temporaldist = (Σ.T .- (ts[indices.-1] - Σ.T * floor.(ts[indices.-1] / Σ.T)))
        @show extrema(temporaldist)
        p = temporaldist ./ dt
        @show extrema(dt)
        dx = xs[:, indices] - xs[:, indices.-1]
        @show extrema(p)
        X = xs[:, indices.-1] + dx * diagm(p)
        t = floor.(ts[indices] / Σ.T) * Σ.T
    end
    if returnindex
        return X, t, indices
    end
    return X, t
end

function section(Σ::SpacialSection, xs, ts; interpolate::Bool=false, returnindex::Bool=false)
    i_before = findall(x -> insideoutside(Σ, x[1], x[2]), collect(zip(eachcol(xs[:, 1:end-2]), eachcol(xs[:, 2:end-1]))))
    i_after = i_before .+ 1
    X = xs[:, i_after]
    t = ts[i_after]

    if interpolate
        H_n = [insideoutside(Σ.plane, col) for col ∈ eachcol(xs[:, i_before])]
        dx_n = [insideoutside(Σ.plane, col) for col ∈ eachcol(xs[:, i_after])] - H_n
        @show extrema(H_n)
        @show extrema(dx_n)
        p = -H_n ./ dx_n
        @show extrema(p)
        dx = xs[:, i_after] - xs[:, i_before]
        X = xs[:, i_before] + dx * diagm(p)
        dt = ts[i_after] - ts[i_before]
        t = ts[i_before] + dt .* p
    end

    if returnindex
        return X, t, i_after
    end
    return X, t
end