# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
using Metis
"""
    part_metis(A, k, alg)

Compute the `k` partions of graph `A` using METIS.

The partitioning method `alg` can be :KWAY or :RECURSIVE.

# Examples
```julia-repl
julia> part_metis(A, 2, :RECURSIVE)
 1
 ⋮
 2
```
"""
function part_metis(A::AbstractSparseMatrix, k::Int, alg::Symbol)

    if alg != :KWAY && alg != :RECURSIVE
        throw(ArgumentError("Invalid algorithm: $alg. Must be :KWAY or :RECURSIVE"))
    end
    
    g = Graph(A)
    p = Metis.partition(g, k, alg = alg)
    p = Int.(p)
    return p
end
