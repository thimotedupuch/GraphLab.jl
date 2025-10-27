# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch


"""
    _recursive_bisection(method::Function, levels::Int, A::AbstractSparseMatrix,
                         coords::Union{Matrix, Nothing}=nothing, minpoints::Int=8,
                         vn::Vector{Int}=Int[])

Recursively partition the graph `A` using the given partitioning `method`, applying hierarchical bisection.

# Arguments
- `method::Function`: Partitioning method to apply (e.g., `part_spectral`, `part_inertial`).
- `levels::Int`: Number of recursive partitioning levels.
- `A::AbstractSparseMatrix`: Adjacency matrix of the graph.
- `coords::Union{Matrix, Nothing}=nothing`: Node coordinates for spatial partitioning (if applicable).
- `minpoints::Int=8`: Minimum number of nodes required to continue partitioning.
- `vn::Vector{Int}=Int[]`: Vector of node indices, used for tracking original node ordering.

# Returns
- A vector of partition labels for each node, recursively refined through hierarchical bisection.

# Example
```julia-repl
julia> _recursive_bisection(part_spectral, 3, A, coords)
 1
 ⋮
 4
 ```
"""
function _recursive_bisection(
    method::Function,
    levels::Int,
    A::AbstractSparseMatrix,
    coords::Union{Matrix,Nothing}=nothing,
    minpoints::Int=8, vn::Vector{Int}=Int[]
)

    n = size(A)[1]

    if isempty(vn)
        vn = collect(1:n)
    end

    if n < minpoints || levels < 1
        return ones(Int, n)
    else
        if !isnothing(coords)
            p = method(A, coords)
            idx1 = findall(x -> x == 1, p)
            idx2 = findall(x -> x == 2, p)
            coords1 = coords[idx1, :]
            coords2 = coords[idx2, :]
        else
            p = method(A)
            idx1 = findall(x -> x == 1, p)
            idx2 = findall(x -> x == 2, p)
            coords1 = coords2 = nothing
        end

        vn1 = vn[idx1]
        vn2 = vn[idx2]

        A1 = A[idx1, idx1]
        A2 = A[idx2, idx2]

        # if !isemtpy(coords)
        if levels > 1
            levels = levels - 1
            p1 = _recursive_bisection(method, levels, A1, coords1, minpoints, vn1)
            p2 = _recursive_bisection(method, levels, A2, coords2, minpoints, vn2)

            return vcat(p1, p2 .+ maximum(p1))[sortperm(vcat(vn1, vn2))]

        end

        return p[sortperm(vn)]
        # end

    end
end


"""
    recursive_bisection(method::Function, k::Int, A::AbstractSparseMatrix,
                        coords::Union{Matrix, Nothing}=nothing, minpoints::Int=8)

Partition the graph `A` into `k` parts using recursive bisection with the specified partitioning `method`.

# Arguments
- `method::Function`: Partitioning method to apply at each bisection step (e.g., `part_spectral`, `part_inertial`).
- `k::Int`: Number of partitions (must be a power of 2 or will be rounded up).
- `A::AbstractSparseMatrix`: Adjacency matrix of the graph.
- `coords::Union{Matrix, Nothing}=nothing`: Node coordinates for spatial partitioning (optional).
- `minpoints::Int=8`: Minimum number of nodes required for further partitioning.

# Returns
- A vector of partition labels for each node.

# Example
```julia-repl
julia> recursive_bisection(part_spectral, 4, A, coords)
 1
 ⋮
 4
```
"""
function recursive_bisection(
    method::Function,
    k::Int,
    A::AbstractSparseMatrix,
    coords::Union{Matrix,Nothing}=nothing,
    minpoints::Int=8
)

    levels = log2(k)
    if !isinteger(levels)
        @warn "log2($k) is not an integer. Rounding up to the closest integer."
    end

    levels = ceil(Int, levels)

    # we will have 2^levels number of partitions.
    return _recursive_bisection(method, levels, A, coords, minpoints)
end
