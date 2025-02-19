module GraphPartitioning

# Load required libraries
using LinearAlgebra
using SparseArrays
using Statistics

# Include all tools
include("tools/recursive_bisection.jl")
include("tools/part_spectral.jl")
include("tools/part_nonspectral.jl")
include("tools/graph_io.jl")
include("tools/part_metis.jl")
include("tools/part_kahip.jl")
include("tools/part_graclus.jl")
include("tools/part_mdc.jl")
include("tools/draw.jl")
include("tools/util.jl")

end # module GraphPartitioning

