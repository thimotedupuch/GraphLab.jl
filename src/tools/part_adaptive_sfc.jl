# Initialize a mutable counter (as a reference)
const traversal_counter = Ref(1);


"""
    KdNode

A node in a 2D kd-tree data structure, representing a spatial region and its recursive subdivision.

Each node contains bounding box coordinates, optional splitting information, and links to child nodes. The node can also store associated data points and optional entry/exit traversal symbols for space-filling curve traversal.

# Fields
- `bbox::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}`: Bounding box as `(min_coords, max_coords)`.
- `splitdim::Union{Int, Nothing}`: Splitting dimension (1 for x-axis, 2 for y-axis), or `nothing` for leaf.
- `splitval::Union{Float64, Nothing}`: Coordinate value of the splitting plane, or `nothing` for leaf.
- `left::Union{KdNode, Nothing}`: Left child node (subtree with coordinates ≤ `splitval`), or `nothing`.
- `right::Union{KdNode, Nothing}`: Right child node (subtree with coordinates > `splitval`), or `nothing`.
- `is_leaf::Bool`: Whether the node is a leaf.
- `is_root::Bool`: Whether the node is the root of the tree.
- `data::Matrix{Float64}`: Data points contained in this node.
- `entry::Union{Symbol, Nothing}`: Optional entry symbol for space-filling curve traversal.
- `exit::Union{Symbol, Nothing}`: Optional exit symbol for space-filling curve traversal.
"""
mutable struct KdNode
    bbox::Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}  # (min_coords, max_coords)
    splitdim::Union{Int, Nothing}
    splitval::Union{Float64, Nothing}
    left::Union{KdNode, Nothing}
    right::Union{KdNode, Nothing}
    is_leaf::Bool
    is_root::Bool
    data::Matrix{Float64}
    entry::Union{Symbol, Nothing}
    exit::Union{Symbol, Nothing}
end


"""
    part_adaptive_sfc(A::SparseMatrixCSC, coords::Matrix, k::Int=2)

Partition a graph using an adaptive space-filling curve traversal.

This function builds a kd-tree from node coordinates, traverses the tree following
an adaptive space-filling curve, and assigns partition labels based on traversal order.
It returns a partition vector that assigns each node to one of `k` partitions.

This algorithm is based on "A General Space-filling Curve Algorithm for Partitioning 2D Meshes" (Aparna et al., 2015, DOI: 10.1109/HPCC-CSS-ICESS.2015.192)
and "Space-filling Curves for Partitioning Adaptively Refined Meshes" (Sasidharan, Aparna, and Snir).


# Arguments
- `A::SparseMatrixCSC`: Adjacency matrix of the graph.
- `coords::Matrix`: An `n × 2` matrix of node coordinates.
- `k::Int=2`: Number of partitions.

# Returns
- A vector `part` of length `n` assigning each node to a partition labeled `1` to `k`.

# Example
```julia-repl
julia> part = part_adaptive_sfc(A, coords, 4)
[1, 1, 2, 2, 3, 4, 4, ...]
"""
function part_adaptive_sfc(A::SparseMatrixCSC, coords::Matrix, k::Int=2)
    n = size(A, 1);

    coords_with_idx = hcat(coords, collect(1:n)); # Assuming that nodes ids are {1, 2, ..., n}

    # Build the kd-Tree
    tree = build_tree(coords_with_idx);

    # Start traversal of the tree
    traversal_counter[] = 1;
    traversal_data = traverse_kd(tree, :L, :R);

    # Extract x, y coordinates and traversal order from each leaf’s data
    xs = [leaf_data[1] for leaf_data in traversal_data];
    ys = [leaf_data[2] for leaf_data in traversal_data];
    sfc_order = Int.([leaf_data[3] for leaf_data in traversal_data]);

    part = get_partition(sfc_order, k);

    return part;
end


"""
    bounding_box(coords::Matrix{Float64})

Compute the axis-aligned bounding box for a set of 2D points.

# Arguments
- `coords::Matrix{Float64}`: A matrix of size `n × 2`, where each row represents a 2D point `(x, y)`.

# Returns
- A tuple `((xmin, ymin), (xmax, ymax))` representing the lower-left and upper-right corners of the bounding box.

# Example
```julia-repl
julia> bounding_box([1.0 2.0; 3.0 4.0; -1.0 5.0])
((-1.0, 2.0), (3.0, 5.0))
"""
function bounding_box(coords::Matrix{Float64})
    xmin = minimum(coords[:, 1])
    ymin = minimum(coords[:, 2])
    xmax = maximum(coords[:, 1])
    ymax = maximum(coords[:, 2])
    return ((xmin, ymin), (xmax, ymax))
end


"""
    split(coords::Matrix{Float64}, dim::Int)

Split a set of 2D points along a specified dimension at the median.

# Arguments
- `coords::Matrix{Float64}`: A matrix of size `n × 2`, where each row is a 2D point `(x, y)`.
- `dim::Int`: Splitting dimension (`1` for x-axis, `2` for y-axis).

# Returns
- `left::Matrix{Float64}`: Points in the left subset (≤ split value).
- `right::Matrix{Float64}`: Points in the right subset (> split value).
- `splitval::Float64`: The coordinate value at which the split occurs.

# Example
```julia-repl
julia> left, right, splitval = split([1.0 2.0; 3.0 4.0; 0.0 5.0], 1)
([0.0 5.0; 1.0 2.0], [3.0 4.0], 1.0)
"""
function split(coords::Matrix{Float64}, dim::Int)
    perm = sortperm(coords[:, dim])
    sorted = coords[perm, :]
    n = size(coords, 1)
    m = div(n, 2)
    splitval = sorted[m, dim]
    left = sorted[1:m, :]
    right = sorted[m+1:end, :]
    return left, right, splitval
end


"""
    build_tree(coords::Matrix{Float64}, is_root::Bool=true)

Recursively build a 2D kd-tree from a set of points. Leaf nodes contain the data points inside their region.

# Arguments
- `coords::Matrix{Float64}`: A matrix of size `n × 2`, where each row is a 2D point `(x, y)`.
- `is_root::Bool=true`: Whether this node is the root of the tree (automatically set during recursion).

# Returns
- A `KdNode` representing the root of the kd-tree.

# Example
```julia-repl
julia> tree = build_tree([1.0 2.0; 3.0 4.0; 0.0 5.0])
KdNode(...)
"""
function build_tree(coords::Matrix{Float64}, is_root::Bool=true)
    n = size(coords, 1)
    if n <= 1
        bbox = bounding_box(coords)
        return KdNode(bbox, nothing, nothing, nothing, nothing, true, is_root, coords, nothing, nothing)
    end

    bbox = bounding_box(coords) 
    spans = bbox[2] .- bbox[1] 
    splitdim = argmax(spans) 
    left_coords, right_coords, splitval = split(coords, splitdim)
    
    left_node = build_tree(left_coords, false)
    right_node = build_tree(right_coords, false)

    return KdNode(bbox, splitdim, splitval, left_node, right_node, false, is_root, zeros(0,3), nothing, nothing)
end


"""
    traverse_kd(node::KdNode, entry::Symbol=:L, exit::Symbol=:R, last_exits::Vector{Any}=[], acc::Vector{Any}=[])

Traverse a kd-tree following an adaptive space-filling curve and collect leaf data.

This function recursively traverses a kd-tree according to entry and exit directions,
determining the order of child traversal at each node using adaptive rules.
At each leaf node, it collects the stored data into an accumulator.

# Arguments
- `node::KdNode`: The kd-tree node to traverse.
- `entry::Symbol=:L`: Entry direction (`:L`, `:R`, `:T`, `:B`) at the current node.
- `exit::Symbol=:R`: Exit direction at the current node.
- `last_exits::Vector{Any}=[]`: List of exit points accumulated during traversal.
- `acc::Vector{Any}=[]`: Accumulator for collected leaf data.

# Returns
- If `node.is_root`, returns `acc`, the accumulated leaf data in traversal order.
- Otherwise returns the last exit point from the traversal.

# Example
```julia-repl
julia> data = traverse_kd(tree)
[[x1, y1, id1], [x2, y2, id2], ...]
"""
function traverse_kd(node::KdNode, entry::Symbol=:L, exit::Symbol=:R, last_exits::Vector{Any}=[], acc::Vector{Any} = [])
    node.entry = entry
    node.exit = exit

    # If the node is empty, stop.
    # If it’s a leaf, count it, save its data in the accumulator, and return its coordinates.
    if node === nothing
        return
    elseif node.is_leaf
        exit_pt = node.bbox[1]
        traversal_counter[] += 1
        push!(acc, node.data)
        return exit_pt
    end

    # Use the node center as entry point if no previous exits.
    # Otherwise reuse the last exit point
    if isempty(last_exits)
        entry_pt = bbox_center(node.bbox)
    else
        entry_pt = last(last_exits)
    end

    # Determine the order of child traversal
    order = step(node.entry, node.exit, node.splitdim, entry_pt, node.splitval)

    # Traverse the first (defined as left) child: 
    # determine its entry/exit directions, recurse, and store its exit point.
    child_l = get_child(node, first(order))
    child_l_entry, child_l_exit = get_child_entry_exit(node.entry, node.exit, node.splitdim, first(order))
    last_exit = traverse_kd(child_l, child_l_entry, child_l_exit, last_exits, acc)
    if !isnothing(last_exit) 
        push!(last_exits, last_exit)
    end
    # Then traverse the second (defined as right) child: 
    # infer its entry as opposite of first child's exit, recurse, and store its exit point.
    child_r = get_child(node, last(order))
    child_r_entry, child_r_exit = (opposite_dir(child_l_exit), node.exit)
    last_exit = traverse_kd(child_r, child_r_entry, child_r_exit, last_exits, acc)
    if !isnothing(last_exit) 
        push!(last_exits, last_exit)
    end

    # Return accumulated data if this is the root node
    if node.is_root
        return acc
    end
end


"""
    bbox_center(bbox)

Compute the center point of a bounding box.

# Arguments
- `bbox`: A tuple `((xmin, ymin), (xmax, ymax))` representing the bounding box corners.

# Returns
- A tuple `(xcenter, ycenter)` representing the center coordinates of the bounding box.

# Example
```julia-repl
julia> bbox_center(((0.0, 0.0), (2.0, 4.0)))
(1.0, 2.0)
"""
function bbox_center(bbox)
    (xmin, ymin), (xmax, ymax) = bbox
    return ((xmin + xmax) / 2, (ymin + ymax) / 2)
end


"""
    step(entry::Symbol, exit::Symbol, splitdim::Int, entry_pt::Tuple{Float64, Float64}, splitval::Float64)

Determine the traversal order of child nodes in a kd-tree based on entry and exit directions.

# Arguments
- `entry::Symbol`: Entry direction (`:L`, `:R`, `:T`, or `:B`).
- `exit::Symbol`: Exit direction (`:L`, `:R`, `:T`, or `:B`).
- `splitdim::Int`: Splitting dimension (`1` for x-axis, `2` for y-axis).
- `entry_pt::Tuple{Float64, Float64}`: Coordinates of the entry point.
- `splitval::Float64`: Value of the splitting plane along `splitdim`.

# Returns
- A vector of symbols indicating the traversal order of child nodes (e.g., `[:left, :right]`).

# Example
```julia-repl
julia> step(:L, :R, 1, (0.0, 0.5), 0.3)
[:left, :right]
"""
function step(entry::Symbol, exit::Symbol, splitdim::Int, entry_pt::Tuple{Float64, Float64}, splitval::Float64)
    type = order_type(entry, exit)
    align = get_alignment(entry, splitdim)

    # Resolve child order based on entry point and split
    if type == :cis
        if align == :parallel
            if splitdim == 1
                return entry == :L ? [:left, :right] : [:right, :left]
            else
                return entry == :T ? [:top, :bottom] : [:bottom, :top]
            end
        else # :cis + :perpendicular
            if splitdim == 1
                return exit == :L ? [:right, :left] : [:left, :right]
                
            else
                return exit == :T ? [:bottom, :top] : [:top, :bottom]
            end
        end
    else  # :tran
        if align == :parallel
            if splitdim == 1 
                return entry == :L ? [:left, :right] : [:right, :left]
            else 
                return entry == :T ? [:top, :bottom] : [:bottom, :top]
            end
        else # :tran + :perpendicular
            if splitdim == 1
                return entry_pt[splitdim] ≤ splitval ? [:left, :right] : [:right, :left]
            else
                return entry_pt[splitdim] ≤ splitval ? [:bottom, :top] : [:top, :bottom]
            end
        end
    end
end


"""
    get_child(node::KdNode, side::Symbol)

Return the child node of a `KdNode` corresponding to the specified side.

This function retrieves either the left or right child node depending on the side and the splitting dimension.

# Arguments
- `node::KdNode`: The kd-tree node.
- `side::Symbol`: Side indicator (`:left`, `:right`, `:bottom`, or `:top`).

# Returns
- The child `KdNode` corresponding to the requested side.

# Example
```julia-repl
julia> get_child(node, :left)
KdNode(...)
"""
function get_child(node::KdNode, side::Symbol)
    if node.splitdim == 1
        return side == :left ? node.left : node.right
    elseif node.splitdim == 2
        return side == :bottom ? node.left : node.right
    else
        error("Invalid splitdim")
    end
end


"""
    get_child_entry_exit(parent_entry::Symbol, parent_exit::Symbol, splitdim::Int, side::Symbol)

Compute the entry and exit directions for a child node in a kd-tree traversal.

Given the parent node’s entry and exit directions, the splitting dimension, and the side being traversed,
this function determines the entry and exit directions for the child node according to space-filling
curve traversal rules (using `:cis` and `:tran` order types and `:parallel` / `:perpendicular` alignment).

# Arguments
- `parent_entry::Symbol`: Entry direction at the parent node (`:L`, `:R`, `:T`, `:B`).
- `parent_exit::Symbol`: Exit direction at the parent node (`:L`, `:R`, `:T`, `:B`).
- `splitdim::Int`: Splitting dimension (`1` for x-axis, `2` for y-axis).
- `side::Symbol`: Side of the child being traversed (`:left`, `:right`, `:top`, `:bottom`).

# Returns
- A tuple `(child_entry, child_exit)` indicating the entry and exit directions for the child node.

# Example
```julia-repl
julia> get_child_entry_exit(:L, :R, 1, :left)
(:L, :R)
"""
function get_child_entry_exit(parent_entry::Symbol, parent_exit::Symbol, splitdim::Int, side::Symbol)
    type = order_type(parent_entry, parent_exit)
    align = get_alignment(parent_entry, splitdim)

    if type == :tran && align == :parallel
        return (parent_entry, parent_exit)
    elseif  type == :tran && align == :perpendicular
        if splitdim == 1
            if parent_entry == :T
                if side == :left
                    return (:T, :R)
                else
                    return (:T, :L)
                end
            else
                if side == :left
                    return (:B, :R)
                else
                    return (:B, :L)
                end
            end
        else
            if parent_entry == :L
                if side == :top
                    return (:L, :B)
                else
                    return (:L, :T)
                end
            else
                if side == :top
                    return (:R, :B)
                else
                    return (:R, :T)
                end
            end

        end
    elseif type == :cis && align == :parallel
        if splitdim == 1
            if side == :left
                if parent_entry == :L
                    return (:L, :R)
                else
                    return (:R, parent_exit)
                end
            else
                if parent_entry == :L
                    return (:L, parent_exit)
                else
                    return (:R, :L)
                end
            end
        else
            if side == :top
                if parent_entry == :T
                    return (:T, :B)
                else
                    return (:B, parent_exit)
                end
            else
                if parent_entry == :T
                    return (:T, parent_exit)
                else
                    return (:B, :T)
                end
            end
        end
    elseif type == :cis && align == :perpendicular
        if splitdim == 1
            if side == :left
                if parent_entry == :T
                    return (:T, :R)
                else
                    return (:B, :R)
                end
            else
                if parent_entry == :T
                    return (:T, :L)
                else
                    return (:B, :L)
                end
            end
        else
            if side == :top
                if parent_entry == :L
                    return (:L, :B)
                else
                    return (:R, :B)
                end
            else
                if parent_entry == :L
                    return (:L, :T)
                else
                    return (:R, :T)
                end
            end
        end
    end
end


"""
    order_type(entry::Symbol, exit::Symbol)

Determine the traversal type based on entry and exit directions.

# Arguments
- `entry::Symbol`: Entry direction (`:L`, `:R`, `:T`, or `:B`).
- `exit::Symbol`: Exit direction (`:L`, `:R`, `:T`, or `:B`).

# Returns
- `:tran` if entry and exit are opposite directions.
- `:cis` otherwise.

# Example
```julia-repl
julia> order_type(:L, :R)
:tran
julia> order_type(:T, :R)
:cis
"""
function order_type(entry::Symbol, exit::Symbol)
    if (entry == :L && exit == :R) || (entry == :R && exit == :L) ||
        (entry == :T && exit == :B) || (entry == :B && exit == :T)
         return :tran
     else
         return :cis
     end
end


"""
    get_alignment(entry::Symbol, splitdim::Int)

Determine the alignment of an entry direction relative to the splitting dimension.

# Arguments
- `entry::Symbol`: Entry direction (`:L`, `:R`, `:T`, or `:B`).
- `splitdim::Int`: Splitting dimension (`1` for x-axis, `2` for y-axis).

# Returns
- `:parallel` if entry edge direction is parallel to the splitting dimension.
- `:perpendicular` otherwise.

# Example
```julia-repl
julia> get_alignment(:L, 1)
:parallel
julia> get_alignment(:T, 1)
:perpendicular
"""
function get_alignment(entry::Symbol, splitdim::Int)
    if (splitdim == 1 && (entry == :L || entry == :R)) ||
       (splitdim == 2 && (entry == :T || entry == :B))
        return :parallel
    else
        return :perpendicular
    end
end


"""
    opposite_dir(dir::Symbol)::Symbol

Return the opposite direction symbol.

# Arguments
- `dir::Symbol`: Direction symbol (`:L`, `:R`, `:T`, or `:B`).

# Returns
- The opposite direction symbol.

# Example
```julia-repl
julia> opposite_dir(:L)
:R
julia> opposite_dir(:T)
:B
"""
function opposite_dir(dir::Symbol)::Symbol
    if dir == :T
        return :B
    elseif dir == :B
        return :T
    elseif dir == :L
        return :R
    elseif dir == :R
        return :L
    else
        error("Unknown direction: $dir")
    end
end


"""
    get_partition(v::AbstractVector, k::Int)

Create a partition indicator vector by dividing indices into `k` contiguous groups.

This function assigns partition labels `1` to `k` to the indices in `v`, distributing them
evenly into `k` partitions based on their position in `v`. Each value in `v` is assumed to
be a node ID (1-based indexing), and the result maps each node ID to its partition.

# Arguments
- `v::AbstractVector`: A vector of node IDs to partition.
- `k::Int`: The number of partitions.

# Returns
- A vector `part` of length `maximum(v)` where `part[id]` is the partition label (1 to `k`) assigned to node `id`.

# Example
```julia-repl
julia> get_partition([1, 2, 3, 4, 5, 6], 2)
[1, 1, 1, 2, 2, 2]
"""
function get_partition(v::AbstractVector, k::Int)
    n = length(v)
    part = zeros(Int, maximum(v))  # assuming v contains node IDs (1-based)

    chunk_size = ceil(Int, n / k)

    for i in 1:k
        start_idx = (i - 1) * chunk_size + 1
        end_idx = min(i * chunk_size, n)
        for j in start_idx:end_idx
            part[v[j]] = i
        end
    end

    return part
end
