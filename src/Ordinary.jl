using Revise

module Ordinary

include("./utils.jl")

include("./Solvers.jl")
include("./Problems.jl")
include("./Sections.jl")
include("./Differences.jl")
include("./DelayedEmbedding.jl")
include("./Lyapunov.jl")
include("./Dimension.jl")

if !isa(CosineDriver(1, 1), AbstractForcing)
    println("ERROR: !isa(CosineDriver(1, 1), AbstractForcing)")
end

end # module Ordinary
