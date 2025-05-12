# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 


"""
    part_coordinate(A::SparseMatrixCSC, coords::Matrix)

Compute a bi-partition of the graph `A` using the coordinate method based on `coords`.

# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix of the graph.
- `coords::Matrix`: Node coordinates used for partitioning.

# Returns
- A vector of partition labels for each node.

# Example
```julia-repl
julia> part_coordinate(A, coords)
 1
 ⋮
 2
```
"""
function  part_coordinate(A::SparseMatrixCSC, coords::Matrix)
    d = size(coords)[2]
    best_cut = Inf
    p = ones(Int, size(coords)[1])
    p_temp = ones(Int, size(coords)[1])
    for dim in 1:d
        v = zeros(d)
        v[dim] = 1
        p1, p2 = _partition(coords, v)
        p_temp[p1] .= 1
        p_temp[p2] .= 2
        this_cut = count_edge_cut(A, p_temp)
        if this_cut < best_cut
            best_cut = this_cut
            p = ones(Int, size(coords)[1])
            p[p1] .= 1
            p[p2] .= 2
        end
    end
    
    return p
end