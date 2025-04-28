# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
using Metis
"""
    part_metis(A::SparseMatrixCSC, k::Int, alg::Symbol)

Partition the graph `A` into `k` parts using METIS with the specified algorithm.

# Arguments
- `A`: Adjacency matrix of the graph.
- `k`: Number of partitions.
- `alg`: Partitioning algorithm (`:KWAY` or `:RECURSIVE`).

# Output
- Returns a vector of partition labels for each node.

# Examples
```julia-repl
julia> part_metis(A, 2, :RECURSIVE)
 1
 ⋮
 2
```
"""
function part_metis(A::SparseMatrixCSC, k::Int=2, alg::Symbol=:KWAY)

    if alg != :KWAY && alg != :RECURSIVE
        throw(ArgumentError("Invalid algorithm: $alg. Must be :KWAY or :RECURSIVE"))
    end
    
    g = Graph(A)
    p = Metis.partition(g, k, alg = alg)
    p = Int.(p)
    return p
end
