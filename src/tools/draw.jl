# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 

using CairoMakie
using Colors

## Note:
# This code is adapted from https://github.com/fcdimitr/SGtSNEpi.jl/blob/master/src/vis.jl
# Here, we augment it for our specific purposes.
# As required by the authors of the original paper, please cite:
# @inproceedings{pitsianis2019sgtsnepi,
#     author = {Pitsianis, Nikos and Iliopoulos, Alexandros-Stavros and Floros, Dimitris and Sun, Xiaobai},
#     doi = {10.1109/HPEC.2019.8916505},
#     booktitle = {IEEE High Performance Extreme Computing Conference},
#     month = {11},
#     title = {{Spaceland Embedding of Sparse Stochastic Graphs}},
#     year = {2019}
# }
"""
    _vis_graph(A::SparseMatrixCSC, coords::Matrix, p::Vector{Int})

Visualize the partitioning `p` of graph `A` using node coordinates `coords`.

# Arguments
- `A`: Adjacency matrix of the graph.
- `coords`: Node coordinates for plotting.
- `p`: Partition labels for each node.

# Output
- Displays the partitioned graph visualization.
"""
function _vis_graph(A::SparseMatrixCSC,
                    coords::Matrix,
                    p::Vector{Int})

    cmap = distinguishable_colors(
      maximum(p) - minimum(p) + 1, # Adjust colormap to the range of labels
      [RGB(1,1,1), RGB(0,0,0)],
      dropseed = true
    )
    
    intra_edge_width = 1.0  # Line width for intra-cluster edges
    inter_edge_width = 0.2  # Line width for inter-cluster edges
    edge_opacity = 1.0      # Transparency for edges
    inter_edge_color = Colors.colorant"gray" # Color for inter-cluster edges
    intra_edge_color = nothing        # Color for intra-cluster edges
    marker_size = 1              # Marker size for scatter points

    # Function to plot lines (edges) for given indices
    function _plot_lines!(ax, Y, i, j, idx, color, lwd)
        ii = vec([reshape(Y[vcat(i[idx], j[idx]), 1], (Int32.(sum(idx)), 2)) NaN * zeros(sum(idx))]')
        jj = vec([reshape(Y[vcat(i[idx], j[idx]), 2], (Int32.(sum(idx)), 2)) NaN * zeros(sum(idx))]')
        lines!(ax, ii, jj, color = color, linewidth = lwd)
    end
        
    # Compute min and max values for x and y
    x_min, x_max = extrema(coords[:, 1])  # x-values
    y_min, y_max = extrema(coords[:, 2])  # y-values

    # Compute aspect ratio
    margin_ratio = 0.05
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_min -= margin_ratio * x_range
    x_max += margin_ratio * x_range
    y_min -= margin_ratio * y_range
    y_max += margin_ratio * y_range
    aspect_ratio = x_range / y_range
    
    height = 800
    width = round(Int, aspect_ratio * height)

    n = size(coords, 1)  # Number of points

    # Normalize labels
    p = p .- minimum(p)
    L_u = sort(unique(p))  # Unique labels

    # Create a figure with specified resolution
    f = Figure(size = (width, height), backgroundcolor=:transparent)

    # Set a global theme for all axes
    set_theme!(Axis = (xgridvisible = false, ygridvisible = false))

    # Create axis
    ax = Axis(f[1, 1],backgroundcolor=:transparent)

    # Plot edges
    i, j = findnz(tril(A))  # Find non-zero elements in the lower triangle of A
    for k in L_u
        idx_inner = map((x, y) -> x == y && x == k, p[i], p[j])
        if isnothing(intra_edge_color)  # Use colormap for intra-cluster edges
            _plot_lines!(ax, coords, i, j, idx_inner, RGBA(cmap[k + 1], edge_opacity), intra_edge_width)
        else
            _plot_lines!(ax, coords, i, j, idx_inner, intra_edge_color, intra_edge_width)
        end
    end
    idx_cross = map((x, y) -> x != y, p[i], p[j])  # Inter-cluster edges
    _plot_lines!(ax, coords, i, j, idx_cross, inter_edge_color, inter_edge_width)

    # Plot scatter points
    scatter!(ax, coords[:, 1], coords[:, 2], color = p, colormap = cmap, markersize = marker_size)
    #ax.aspect = DataAspect()  # Keep aspect ratio consistent
    
    # tightlimits!(ax)          # Fit the view to the data
    limits!(ax, x_min, x_max, y_min, y_max)
    resize_to_layout!(f)
    #ax.aspect = DataAspect()
    return f
end


"""
    draw_graph(A::SparseMatrixCSC, coords::Matrix, p::Vector{Int}; file_name::Union{String, Nothing}=nothing)

Draw and optionally save a visualization of the partitioned graph `A`.

# Arguments
- `A`: Adjacency matrix of the graph.
- `coords`: Node coordinates for plotting.
- `p`: Partition labels for each node.
- `file_name` (optional): If provided, saves the figure to the specified file.

# Output
- Returns the generated figure.
- Saves the figure if `file_name` is specified.
"""
function draw_graph(A::SparseMatrixCSC, coords::Matrix, p::Vector{Int};file_name::Union{String, Nothing}=nothing)
    CairoMakie.activate!()
    fig = _vis_graph(A, coords, p)

    # Save the figure if file_name is provided
    if file_name !== nothing
        save(file_name, fig)
        println("Figure saved as '$file_name'")
    end

    return fig
end


"""
    draw_graph(A::SparseMatrixCSC, coords::Matrix; file_name::Union{String, Nothing}=nothing)

Draw and optionally save a visualization of the graph `A` without partitioning.

# Arguments
- `A`: Adjacency matrix of the graph.
- `coords`: Node coordinates for plotting.
- `file_name` (optional): If provided, saves the figure to the specified file.

# Output
- Returns the generated figure.
- Saves the figure if `file_name` is specified.
"""
function draw_graph(A::SparseMatrixCSC, coords::Matrix; file_name::Union{String, Nothing}=nothing)
    n = size(A)[1]
    p = ones(Int, n)
    return draw_graph(A, coords, p, file_name=file_name)
end


"""
    spy_mat(A::SparseMatrixCSC)

Generates a visual sparsity plot of the sparse matrix `A`.

# Arguments
- `A::SparseMatrixCSC`: The sparse matrix to visualize.

# Returns
- `f::Figure`: The generated `Figure` object containing the sparsity plot.

# Notes
- Uses reversed y-axis to match standard matrix visualization conventions.
"""

function spy_mat(A::SparseMatrixCSC)
    f = Figure(resolution = (1000, 1000), backgroundcolor = :transparent)

    ax1 = Axis(f[1, 1],
        aspect = 1,
        title = "Original\nnnz = $(nnz(A))",
        titlesize = 14,
        titlecolor = :white,
        backgroundcolor = :transparent,
        yreversed = true
    )

   
    # Plot sparsity patterns
    spy!(ax1, A, framecolor = :lightgrey, color = :white)
 

    # Clean up axis decorations
    #foreach(hidedecorations!, (ax1))

    return f
end
