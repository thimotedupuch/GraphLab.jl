using Pkg
Pkg.activate(@__DIR__)

using SparseArrays,DelimitedFiles, MAT, Plots, PrettyTables
using GraphPartitioning

"""
    benchmark(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, coords::Matrix{Int64}, prefix::String)

Run various graph partitioning methods on `A` and save visualizations.

# Arguments
- `A::SparseMatrixCSC{Float64, Int64}`: Adjacency matrix of the graph.
- `coords::Matrix{Int64}`: Node coordinates for visualization.
- `prefix::String`: Filename prefix for saved images.

# Methods Used
- Coordinate-based (`part_coordinate`)
- Inertial (`part_inertial`)
- Spectral (`part_spectral`)
- METIS recursive bisection (`part_metis`)

# Output
- Saves partitioned graph visualizations as PNG files.
"""
function benchmark( A::SparseMatrixCSC, 
                    coords::Matrix, 
                    prefix::String)
    # Partition methods
    p_coord    = GraphPartitioning.part_coordinate(A, coords)
    p_inertial = GraphPartitioning.part_inertial(A, coords)
    p_spectral = GraphPartitioning.part_spectral(A)
    p_metis    = GraphPartitioning.part_metis(A, 2, :RECURSIVE)

    # Save visualizations
    GraphPartitioning.draw_graph(A, coords, p_coord, file_name=prefix * "_coordinate.png")
    GraphPartitioning.draw_graph(A, coords, p_inertial, file_name=prefix * "_inertial.png")
    GraphPartitioning.draw_graph(A, coords, p_spectral, file_name=prefix * "_spectral.png")
    GraphPartitioning.draw_graph(A, coords, p_metis, file_name=prefix * "_metis.png")

    # Compute edge cuts
    results = [
        GraphPartitioning.count_edge_cut(A, p_coord),
        GraphPartitioning.count_edge_cut(A, p_inertial),
        GraphPartitioning.count_edge_cut(A, p_spectral),
        GraphPartitioning.count_edge_cut(A, p_metis)
    ]

    return results
end

# List of .mat files
files = [
    "3elt.mat",
    "airfoil1.mat",
    "barth4.mat",
    "crack.mat",
    "mesh1e1.mat",
    "mesh2e1.mat",
    "mesh3e1.mat",
    "netz4504_dual.mat",
    "stufe.mat",
    "ukerbe1.mat"
]

# Directory containing the files
base_dir = "meshes/"

final_results = Matrix{Any}(undef,length(files),5)

# Process each file
for (i, file_name) in enumerate(files)
    # Open the .mat file
    file_path = joinpath(base_dir, file_name)
    file = MAT.matopen(file_path)
    data = MAT.read(file)

    # Extract adjacency matrix and coordinates
    # note "Swiss_graph.mat" as sliglty different data structure
    if file_name == "Swiss_graph.mat" 
        A = sparse(data["CH_adj"])
        coords = data["CH_coords"]
    else
        A = data["Problem"]["A"]
        coords = data["Problem"]["aux"]["coord"]
    end

    # Symmetrize A
    A = (A + transpose(A))/2

    # Generate a prefix for output files based on the file name
    prefix = "ex1_" * replace(file_name, ".mat" => "")

    # Run the benchmark function
    println("Processing file: $file_path with prefix: $prefix")
    results = benchmark(A, coords, prefix)

    final_results[i,1] = file_name
    final_results[i,2] = results[1]
    final_results[i,3] = results[2]
    final_results[i,4] = results[3]
    final_results[i,5] = results[4]
   
    # Close the .mat file
    MAT.close(file)
end

# Pretty print the final results
header =  ["Mesh Name", "part_coordinate", "part_inertial", "part_spectral", "part_metis"]
pretty_table(final_results, header=header)