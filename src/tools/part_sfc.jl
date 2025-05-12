"""
    part_sfc(A::SparseMatrixCSC, coords::Matrix, k::Int=2, curve::Symbol=:gilbert)

Partition the graph `A` into `k` parts using a space-filling curve on `coords`.

The coordinates are projected to a grid and ordered according to the selected
space-filling curve (`:gilbert` or `:morton`). The resulting 1D ordering
is then split into `k` contiguous segments.

# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix of the graph.
- `coords::Matrix`: Node coordinates used for projection to grid.
- `k::Int`: Number of partitions (default = 2).
- `curve::Symbol`: Curve type, one of `:gilbert` or `:morton`.

# Returns
- A vector of partition labels (1 to `k`) for each node.

# Example
```julia-repl
julia> part_sfc(A, coords, 4, :morton)
 1
 ⋮
 4
```
"""
function part_sfc(A::SparseMatrixCSC, coords::Matrix, k::Int=2, curve::Symbol=:gilbert)
    N = size(coords, 1);

    # 1. Coordinates normalization
    coords_centered = coords .- minimum(coords, dims=1);
    scale = maximum(coords_centered);
    norm_coords = coords_centered ./ scale;

    # 2. Discretization into square grid
    grid_size = ceil(Int, sqrt(N));
    grid_coords = clamp.(floor.(Int, norm_coords .* (grid_size - 1)) .+ 1, 1, grid_size);

    # 3. Fill grid with nodes indices
    grid = [Int[] for _ in 1:grid_size, _ in 1:grid_size];
    for node in 1:N
        x, y = grid_coords[node, :]
        push!(grid[x, y], node)
    end

    # 4. Traverse with space filling curve
    order = _sfc_indices(curve, (grid_size, grid_size));

    # Ordering sanity check
    expected_count = grid_size*grid_size
    @assert length(order) == expected_count "Missing or extra points in order"
    @assert length(unique(order)) == expected_count "Duplicate grid cells in order"

    node_order = reduce(vcat, (grid[I] for I in order))

    # 5. Assign partitions
    part = zeros(Int, N)
    chunk = ceil(Int, N / k)
    for (i, node) in enumerate(node_order)
        part[node] = 1 + div(i - 1, chunk)
    end

    return part

end


function _sfc_indices(curve::Symbol, grid_dims::Tuple{Int, Int}) 
    if curve == :gilbert
        return _gilbert_indices(grid_dims)
    elseif curve == :morton 
        return _morton_indices(grid_dims)
    else
        throw(ArgumentError("Unsupported SFC method: $curve. Use :gilbert or :morton."))
    end
end

# The following function was adapted from code by Simon Byrne's GilbertCurves.jl (https://github.com/CliMA/GilbertCurves.jl/)
# and by Jakub Červený (https://github.com/jakubcerveny/gilbert)


"""
    _gilbert_indices(grid_dims::Tuple{Int, Int}; maj_axis=grid_dims[1] ≥ grid_dims[2] ? 1 : 2)

Generate a list of `CartesianIndex{2}` representing a generalized space-filling traversal
of a grid with dimensions `grid_dims`, using a Gilbert-like curve. Traversal will
favor rows (X axis) if `maj_axis == 1`, or columns (Y axis) if `maj_axis == 2`.

### Arguments
- `grid_dims::Tuple{Int, Int}`: a tuple `(rows, cols)` specifying grid size.
- `maj_axis`: optional; either `1` (row-major preference) or `2` (column-major preference).
   - The curve will traverse more linearly along the `maj_axis`.

### Returns
- A `Vector{CartesianIndex{2}}` representing the traversal path.

### Example
```julia-repl
julia> order = _gilbert_indices((4, 5))
julia> println(order)
[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(2, 2), ...]
```
"""
_gilbert_indices(grid_dims::Tuple{Int, Int}; maj_axis=grid_dims[1] ≥ grid_dims[2] ? 1 : 2) =
_gilbert_order(CartesianIndices(grid_dims); maj_axis)


"""
    _gilbert_order(grid_indices::AbstractMatrix; maj_axis=rows ≥ cols ? 1 : 2)

Return a list of `CartesianIndex` values from a grid, traversed recursively using
Gilbert curve logic. When `maj_axis == 2`, the grid is transposed to preserve logic.
"""
function _gilbert_order(grid_indices::AbstractMatrix; maj_axis=size(grid_indices, 1) ≥ size(grid_indices, 2) ? 1 : 2)
    traversal = CartesianIndex{2}[]
    if maj_axis == 1
        _append_gilbert!(traversal, grid_indices)
    else
        _append_gilbert!(traversal, permutedims(grid_indices, (2, 1)))
    end
    return traversal
end


"""
    _append_gilbert!(path, grid::AbstractMatrix)

Recursively traverse the given matrix `grid` using a generalized Gilbert curve pattern.
Appends the `CartesianIndex` entries to `path` in order.
"""
function _append_gilbert!(path::Vector{CartesianIndex{2}}, grid::AbstractMatrix)
    nrows, ncols = size(grid)

    if nrows == 1 || ncols == 1
        append!(path, vec(grid))  # base case: single row or column
    elseif 2*nrows > 3*ncols
        mid = div(nrows, 2)
        if isodd(mid) && nrows > 2
            mid += 1
        end
        _append_gilbert!(path, grid[1:mid, :])
        _append_gilbert!(path, grid[mid+1:end, :])
    else
        row_mid = div(nrows, 2)
        col_mid = div(ncols, 2)
        if isodd(col_mid) && ncols > 2
            col_mid += 1
        end
        # Top-left block, transposed
        _append_gilbert!(path, permutedims(grid[1:row_mid, 1:col_mid], (2, 1)))
        # Right block
        _append_gilbert!(path, grid[:, col_mid+1:end])
        # Bottom-left block, transposed and reversed
        _append_gilbert!(path, permutedims(grid[end:-1:row_mid+1, col_mid:-1:1], (2, 1)))
    end
end


# The following Morton encoding functions are based on the standard
# bit interleaving technique for computing Morton codes (Z-order curve indices).
# This implementation was adapted from:
# - https://en.wikipedia.org/wiki/Z-order_curve
# - libmorton (https://github.com/Forceflow/libmorton)
"""
    _morton_indices(grid_dims::Tuple{Int, Int})

Return a list of `CartesianIndex{2}` grid coordinates ordered by 2D Morton (Z-order) curve.

### Arguments
- `grid_dims::Tuple{Int, Int}`: the number of rows and columns of the grid, e.g. `(64, 64)`.

### Returns
- A vector of `CartesianIndex{2}` objects representing the traversal order.

### Example
```julia-repl
julia> order = _morton_indices((4, 4))
julia> println(order)
[CartesianIndex(1, 1), CartesianIndex(1, 2), CartesianIndex(2, 1), CartesianIndex(2, 2), ...]
```
"""
function _morton_indices(grid_dims::Tuple{Int, Int})
    rows, cols = grid_dims
    order = CartesianIndex{2}[]

    for i in 1:rows, j in 1:cols
        push!(order, CartesianIndex(i, j))
    end

    # Sort by Morton index (bit interleaving of (i-1, j-1))
    sort!(order, by = ij -> _morton_index(ij[1] - 1, ij[2] - 1))
    return order
end


"""
    _morton_index(x::Int, y::Int)

Interleave the bits of x and y to produce a Morton code (Z-order curve).
"""
function _morton_index(x::Int, y::Int)
    return _part1by1(x) | (_part1by1(y) << 1)
end


"""
    _part1by1(n::Int)

Expand 16-bit integer `n` into 32 bits with zeros between the original bits.
"""
function _part1by1(n::Int)
    n = n & 0x0000ffff
    n = (n | (n << 8))  & 0x00FF00FF
    n = (n | (n << 4))  & 0x0F0F0F0F
    n = (n | (n << 2))  & 0x33333333
    n = (n | (n << 1))  & 0x55555555
    return n
end