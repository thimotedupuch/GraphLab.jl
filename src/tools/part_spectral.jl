# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
using Arpack, Statistics, Random


"""
    part_spectral(A::SparseMatrixCSC; fiedler::Bool=false)

Compute a bi-partition of the graph `A` using the spectral method.

# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix of the graph.
- `fiedler::Bool=false`: If `true`, returns the Fiedler vector instead of partition labels.

# Returns
- A vector of partition labels (1 or 2).
- If `fiedler=true`, returns the Fiedler vector.

# Example
```julia-repl
julia> part_spectral(A)
 1
 ⋮
 2
```
"""
function part_spectral(A::SparseMatrixCSC; fiedler::Bool=false)
    n = size(A)[1]

    if n > 4*10^4
        @warn "graph is large. Computing eigen values may take too long."     
    end

    # 1. Construct the Laplacian matrix.
    D = spzeros(n, n)
    degrees = vec(sum(A, dims=1))

    for i in 1:n
        D[i, i] = degrees[i]
    end
    L = D - A

    # 2. Compute its eigendecomposition.
    # v0 = ones(n)
    local_rng = MersenneTwister(1234)
    v0 = randn(local_rng, n)
    eig_vals, eig_vecs = Arpack.eigs(L; which=:SR, nev=2, v0=v0)

    ev2 = eig_vecs[:, sortperm(eig_vals)[2]]
    # ev3 = eig_vecs[:, sortperm(eig_vals)[3]]

    fiedler_vec = ev2
    # 3. Label the vertices with the entries of the Fiedler vector.
    m = median(fiedler_vec)
    # 4. Partition them around their median value, or 0.
    p = map(x->x .> m, fiedler_vec).+1
    # 5. Return the indicator vector
    return p
end
