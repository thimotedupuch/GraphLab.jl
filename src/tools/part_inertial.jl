# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch  
using Random


"""
    part_inertial(A::SparseMatrixCSC, coords::Matrix)

Compute a bi-partition of the graph `A` using the inertial method based on `coords`.

# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix of the graph.
- `coords::Matrix`: Node coordinates used for partitioning.

# Returns
- A vector of partition labels for each node.

# Example
```julia-repl
julia> part_inertial(A, coords)
 1
 ⋮
 2
```
"""
function part_inertial(A::SparseMatrixCSC, coords::Matrix)
    x = coords[:, 1]
    y = coords[:, 2]

    n = size(x)[1]

    # 1. Compute the center of mass.
    x_bar = sum(x)/n
    y_bar = sum(y)/n

    # 2. Construct the matrix M. (see pdf of the assignment)
    sxx = sum((x .- x_bar).^2)
    syy = sum((y .- y_bar).^2)
    sxy = sum((x .- x_bar).*(y .- y_bar))
    M = [[syy sxy]
        [sxy sxx]]

    # 3. Compute the eigenvector associated with the smallest eigenvalue of M.
    # v0 = ones(size(M)[1])
    local_rng = MersenneTwister(1234)
    v0 = randn(local_rng, size(M)[1])
    eig_vals, eig_vecs = eigs(M; which=:SR, nev=1, v0=v0)

    # 4. Partition the nodes around line L 
    #    (use may use the function partition(coords, eigv))
    i = sortperm(eig_vals)
    eigv = eig_vecs[:,i]
    eigv = [eigv[1]; -eigv[2]]
    eigv = eigv./norm(eigv)
    p1, p2 = _partition(coords, eigv)
    p = ones(Int, n)
    p[p1] .= 1
    p[p2] .= 2

    # 5. Return the indicator vector
    return p

end



