using Graphs

"""
    count_edge_cut(A, p)

Counts the number of edges that cross partitions in a graph.

# Arguments
- `A`: Adjacency matrix of the graph (can be sparse).
- `p`: Partition vector where `p[v]` represents the partition of vertex `v`.

# Returns
The number of edges that connect nodes in different partitions.
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

function build_adjacency(edges::Matrix{Int}, num_nodes::Int)
    adjacency_matrix = spzeros(Int, num_nodes, num_nodes)
    for edge in eachrow(edges)
        adjacency_matrix[edge[1], edge[2]] = 1
        adjacency_matrix[edge[2], edge[1]] = 1  # Ensure symmetry
    end
    return adjacency_matrix
end

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
