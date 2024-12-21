# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
using Arpack 

"""
    part_spectral(A, fiedler=false)

Compute the bi-partions of graph `A` using spectral method.

If `fiedler` is true, return the entries of the fiedler vector.

# Examples
```julia-repl
julia> part_spectral(A)
 1
 ⋮
 2
```
"""

function part_spectral(A)
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
    eig_vals, eig_vecs = Arpack.eigs(L; which=:SR, ncv = 300)

    ev2 = eig_vecs[:, sortperm(eig_vals)[2]]
    ev3 = eig_vecs[:, sortperm(eig_vals)[3]]

    fiedler_vec = ev2
    # 3. Label the vertices with the entries of the Fiedler vector.
    # 4. Partition them around their median value, or 0.
    p = map(x->x .> 0, fiedler_vec).+1
    # 5. Return the indicator vector
    return p
end
