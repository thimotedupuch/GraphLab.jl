# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 
"""
    draw_graph(A, coords, p)

Plot the partitions `p` of graph `A`.
"""

using CairoMakie
using Colors

## Note:
# This code is adapted from https://github.com/fcdimitr/SGtSNEpi.jl/blob/master/src/vis.jl
# Here, we augment it for our specific purposes.
# As required by the authors of the original paper, please cite:
# }
function _vis_graph(A, coords, p)
    
    cmap = distinguishable_colors(
      maximum(p) - minimum(p) + 1, # Adjust colormap to the range of labels
      [RGB(1,1,1), RGB(0,0,0)],
      dropseed = true
    )
    res = (800, 800)        # Resolution of the output figure
    intra_edge_width = 1.0  # Line width for intra-cluster edges
    inter_edge_width = 0.1  # Line width for inter-cluster edges
    edge_opacity = 1.0      # Transparency for edges
    inter_edge_color = Colors.colorant"gray" # Color for inter-cluster edges
    intra_edge_color = nothing        # Color for intra-cluster edges
    marker_size = 4                   # Marker size for scatter points

    # Function to plot lines (edges) for given indices
    function _plot_lines!(ax, Y, i, j, idx, color, lwd)
        ii = vec([reshape(Y[vcat(i[idx], j[idx]), 1], (Int32.(sum(idx)), 2)) NaN * zeros(sum(idx))]')
        jj = vec([reshape(Y[vcat(i[idx], j[idx]), 2], (Int32.(sum(idx)), 2)) NaN * zeros(sum(idx))]')
        lines!(ax, ii, jj, color = color, linewidth = lwd)
    end

    n = size(coords, 1)  # Number of points

    # Normalize labels
    p = p .- minimum(p)
    L_u = sort(unique(p))  # Unique labels

    # Create a figure with specified resolution
    f = Figure(size = res,backgroundcolor=:transparent)

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
    tightlimits!(ax)          # Fit the view to the data
    
    return f
end


function _alpha_colorbuffer(figure)
    scene = figure.scene
    bg = scene.backgroundcolor[]
    scene.backgroundcolor[] = RGBAf(0, 0, 0, 1)
    b1 = copy(colorbuffer(scene))
    scene.backgroundcolor[] = RGBAf(1, 1, 1, 1)
    b2 = colorbuffer(scene)
    scene.backgroundcolor[] = bg
    return map(b1, b2) do b1, b2
        calculate_rgba(b1, b2, bg)
    end
end

function draw_graph(A::AbstractSparseMatrix, coords::Matrix, p::Vector{Int};file_name::Union{String, Nothing}=nothing)
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
    draw_graph(A, coords)

Plot the graph based on its adjacency matrix `A`.
"""
function draw_graph(A::AbstractSparseMatrix, coords::Matrix;file_name::Union{String, Nothing}=nothing)
    n = size(A)[1]
    p = ones(Int, n)
    return draw_graph(A, coords, p,file_name=file_name)
end
