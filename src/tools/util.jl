using Graphs


"""
    count_edge_cut(A::AbstractMatrix, p::AbstractVector)

Count the number of edges that cross partitions in the graph `A`.

# Arguments
- `A::AbstractMatrix`: Adjacency matrix of the graph.
- `p::AbstractVector`: Partition vector where `p[v]` represents the partition of vertex `v`.

# Returns
- The number of edges that connect nodes in different partitions.

# Example
```julia-repl
julia> count_edge_cut(A, p)
15
```
"""
function count_edge_cut(A::AbstractMatrix, p::AbstractVector)
    # Convert adjacency matrix to Graph object
    g = Graphs.SimpleGraph(A)

    # Initialize the edge cut count
    edge_cut = 0

    # Iterate over edges
    for e in Graphs.edges(g)
        u, v = Graphs.src(e), Graphs.dst(e)
        if p[u] != p[v]  # Check if edge crosses partitions
            edge_cut += 1
        end
    end

    return edge_cut
end


"""
    build_adjacency(edges::Matrix{Int}, num_nodes::Int)

Construct the adjacency matrix of a graph from an edge list.

# Arguments
- `edges::Matrix`: Matrix where each row represents an edge `[u, v]`.
- `num_nodes::Int`: Total number of nodes in the graph.

# Returns
- A symmetric sparse adjacency matrix (`SparseMatrixCSC{Int, Int}`).

# Example
```julia-repl
julia> edges = [1 2; 2 3; 3 1]
julia> A = build_adjacency(edges, 3)
```
"""
function build_adjacency(edges::Matrix{Int}, num_nodes::Int)
    adjacency_matrix = spzeros(Int, num_nodes, num_nodes)
    for edge in eachrow(edges)
        adjacency_matrix[edge[1], edge[2]] = 1
        adjacency_matrix[edge[2], edge[1]] = 1  # Ensure symmetry
    end
    return adjacency_matrix
end


"""
    build_adjacency(type::String)

Generate a predefined adjacency matrix and corresponding node coordinates.

# Arguments
- `type::String`: Type of graph to generate.  
  - `"network"`: A predefined network structure.  
  - `"triangles"`: A small triangular mesh structure.

# Returns
- `A::SparseMatrixCSC`: The sparse adjacency matrix of the graph.
- `coords::Matrix`: Node coordinates for visualization.

# Example
```julia-repl
julia> A, coords = build_adjacency("network")
```
"""
function build_adjacency(type::String)
    if type=="network"
        A = [[0. 1. 1. 0. 0. 1. 0. 0. 1. 1.]
            [1. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
            [1. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
            [0. 0. 0. 0. 1. 1. 0. 0. 0. 0.]
            [0. 0. 0. 1. 0. 1. 0. 0. 0. 0.]
            [1. 0. 0. 1. 1. 0. 1. 1. 0. 0.]
            [0. 0. 0. 0. 0. 1. 0. 1. 0. 0.]
            [0. 0. 0. 0. 0. 1. 1. 0. 0. 0.]
            [1. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
            [1. 0. 0. 0. 0. 0. 0. 0. 1. 0.]]
        A = sparse(A)
        x = [3 1 2 4 5 3 2 1 5 4]
        y = [2 1 1 5 5 4 5 5 1 1]
    elseif type=="triangles"
        A = [   [0. 1. 1. 1. 0. 0. ]
                [1. 0. 1. 0. 1. 0. ]
                [1. 1. 0. 0. 0. 1. ]
                [1. 0. 0. 0. 1. 1. ]
                [0. 1. 0. 1. 0. 1. ]
                [0. 0. 1. 1. 1. 0. ]]
        A = sparse(A)
        x = [2 4 4 0 6 6]
        y = [3 2 4 3 0 6]
    end
    
    coords = hcat(vec(x), vec(y))

    return A, coords
end
