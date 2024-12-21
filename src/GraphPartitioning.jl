module GraphPartitioning

# Load required libraries
using LinearAlgebra
using SparseArrays
using Statistics


# Include all tools
include("tools/recursive_bisection.jl")
include("tools/part_spectral.jl")
include("tools/part_nonspectral.jl")
include("tools/part_metis.jl")
include("tools/draw.jl")
include("tools/auxiliary.jl")


function test()
    @info "Running GraphPartitioning tests..."
    include("../test/runtests.jl")  # Path relative to the src directory
    @info "All tests completed."
end


end # module GraphPartitioning

