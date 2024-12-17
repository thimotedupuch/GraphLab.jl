# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
"""
    inertial_part(A, coords)

Compute the bi-partions of graph `A` using inertial method based on the
`coords` of the graph.

# Examples
```julia-repl
julia> inertial_part(A, coords)
 1
 ⋮
 2
```
"""
function inertial_part(A, coords)
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
    eig_vals, eig_vecs = eigs(M; which=:SR, ncv = 300, nev=1)

    # 4. Partition the nodes around line L 
    #    (use may use the function partition(coords, eigv))
    i = sortperm(eig_vals)
    eigv = eig_vecs[:,i]
    eigv = [eigv[1]; -eigv[2]]
    eigv = eigv./norm(eigv)
    p1, p2 = partition(coords, eigv)
    p = ones(Int, n)
    p[p1] .= 1
    p[p2] .= 2

    # 5. Return the indicator vector
    return p

end
