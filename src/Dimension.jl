export dcap

using SparseArrayKit

struct Hypercube
    mins::Vector{Float64}
    sides::Vector{Float64}
end

function Hypercube(X::Vector{Vector})
    Hypercube(reduce(hcat, X))
end

function Hypercube(X::Matrix)
    mins, maxs = unzip(extrema(X, dims=(2,)))
    Hypercube(vec(mins), vec(maxs - mins))
end

# helper functions

# Maps x ∈ H ⊂ ℜⁿ to a point in normalized hypercube xₙ ∈ [0,1]ⁿ
function renormalize2hypercube!(H::Hypercube, x)
    (; mins, sides) = H
    x .-= mins
    x ./= sides
    x .= min.(x, 1 - eps())
end

function renormalize2hypercube(H::Hypercube, x)
    min.((x .- H.mins) ./ H.sides, 1 - eps())
end

# Computes the number of ϵ-cubes needed to cover each side of a bounding box
function cubesPerSide(H::Hypercube, ϵ::Float64)
    ceil.(Int, H.sides / ϵ)
end

# Index of normalized point in boolean array (1-indexed)
function hypercubeIndex(x_norm, cubes_per_side)
    Tuple(floor.(Int, x_norm .* cubes_per_side) .+ 1)
end

Base.@inline function N_coarsen(A₀, scale::Real)
    @time length(Set(ceil.(Int, i ./ scale) for i in A₀))
end

function dcap(X_in::Matrix, ϵ₀::Float64; scale_factor::Real=2.0)
    H = Hypercube(X_in)
    @time X = renormalize2hypercube.(Ref(H), eachcol(X_in))
    @time cubes_per_side = cubesPerSide(H, ϵ₀)
    M = floor(Int, log(scale_factor, minimum(cubes_per_side)))
    ϵ = zeros(Float64, M)
    N = zeros(Int, M)
    # A = SparseArray{Bool,length(cubes_per_side)}(undef, tuple(cubes_per_side...))
    # indices = map(x -> hypercubeIndex(renormalize2hypercube(H, x), cubes_per_side), eachcol(X))
    # for i in Set(hypercubeIndex.(eachcol(X), Ref(cubes_per_side)))
    nonzero_indices = Set(hypercubeIndex(x, cubes_per_side) for x in X)
    ϵ[1] = ϵ₀
    N[1] = length(nonzero_indices)
    for m in 2:M
        ϵ[m] = scale_factor^(m - 1) * ϵ₀
        @time N[m] = N_coarsen(nonzero_indices, scale_factor^(m - 1))
        # if abs(ϵ[m] - 0.5756503766994984) < 1e-6
        #     open("indices.jl", "w") do io
        #         println(io, "ϵ = ", ϵ[m])
        #         println(io)
        #         println(io, "indices = [")
        #         println.(io, Ref("    "), Set(ceil.(Int, i ./ scale_factor^(m - 1)) for i in nonzero_indices), Ref(","))
        #         println(io, "]")
        #     end
        # end
    end
    ϵ, N
end