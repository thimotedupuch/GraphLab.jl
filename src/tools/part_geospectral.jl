# ------------------------------------------------------------------------
# Based on work by:
# - John R. Gilbert (UC Santa Barbara): http://www.cs.ucsb.edu/~gilbert/
# - Shang-Hua Teng (University of Southern California): http://www-bcf.usc.edu/~shanghua/
# - Yingzhou Li (Fudan University): https://yingzhouli.com
#
# Note:
# The following acknowledgement is from Yingzhou Li's original codebase:
# "Thanks to Tim Davis for updating the toolbox to MATLAB 5"

using Arpack


"""
    part_geospectral(A::SparseMatrixCSC; ev::Int=2)

Performs spectral partitioning by projecting the graph onto its low-frequency Laplacian eigenvectors,
followed by geometric random sphere partitioning in the resulting embedding space.

# Arguments
- `A::SparseMatrixCSC`: The adjacency matrix of an undirected graph.
- `ev::Int=2`: The number of Laplacian eigenvectors to use for embedding (default is 2).

# Returns
- `p::Vector{Int}`: A partition of the vertex set, typically as cluster labels.

# Notes
- Emits a warning if the graph is large, as eigen decomposition may become computationally expensive.
- This combines spectral embedding with geometric partitioning for improved structural fidelity.
"""
function part_geospectral(A::SparseMatrixCSC; ev::Int=2)
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
    _, eig_vecs = Arpack.eigs(L; which=:SR, nev=ev, ncv = 300)

    # 3. Applies geometric partitioning in the embedding space.
    p = part_randsphere(A, eig_vecs)

    return p
end