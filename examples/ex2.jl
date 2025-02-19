using Pkg
Pkg.activate(@__DIR__)

using SparseArrays,DelimitedFiles, MAT, Plots, PrettyTables
using GraphPartitioning

"""
    benchmark(A::SparseMatrixCSC{Float64, Int64}, coords::Matrix{Int64}, k::Int64, prefix::String)

Run various graph partitioning methods on `A` and save visualizations.

# Arguments
- `A::SparseMatrixCSC{Float64, Int64}`: Adjacency matrix of the graph.
- `coords::Matrix{Int64}`: Node coordinates for visualization.
- `k::Int64`: Number of partitions.
- `prefix::String`: Filename prefix for saved images.

# Methods Used
- Coordinate-based (`part_coordinate`)
- Inertial (`part_inertial`)
- Spectral (`part_spectral`)
- METIS recursive bisection (`part_metis` with `k` partitions).

# Output
- Saves partitioned graph visualizations as PNG files.
"""
function benchmark( A::SparseMatrixCSC, 
                    coords::Matrix, 
                    k::Int64, 
                    prefix::String)


    # Partition methods
    p_coord    = GraphPartitioning.recursive_bisection(GraphPartitioning.part_coordinate,k,A, coords)
    p_inertial = GraphPartitioning.recursive_bisection(GraphPartitioning.part_inertial,  k,A,coords)
    p_spectral = GraphPartitioning.recursive_bisection(GraphPartitioning.part_spectral,  k,A)

    p_metis_rec     = GraphPartitioning.part_metis(A, k, :RECURSIVE)
    p_metis_kway    = GraphPartitioning.part_metis(A, k, :KWAY)

    # Save visualizations
    GraphPartitioning.draw_graph(A, coords, p_coord, file_name=prefix * "_coordinate.png")
    GraphPartitioning.draw_graph(A, coords, p_inertial, file_name=prefix * "_inertial.png")
    GraphPartitioning.draw_graph(A, coords, p_spectral, file_name=prefix * "_spectral.png")
    GraphPartitioning.draw_graph(A, coords, p_metis_rec, file_name=prefix * "_metis_rec.png")
    GraphPartitioning.draw_graph(A, coords, p_metis_kway, file_name=prefix * "_metis_kway.png")

    # Compute edge cuts
    results = [
        GraphPartitioning.count_edge_cut(A, p_coord),
        GraphPartitioning.count_edge_cut(A, p_inertial),
        GraphPartitioning.count_edge_cut(A, p_spectral),
        GraphPartitioning.count_edge_cut(A, p_metis_rec),
        GraphPartitioning.count_edge_cut(A, p_metis_kway)
    ]

    return results
end

# List of .mat files
files = [
    "Swiss_graph.mat"
]

# Directory containing the files
base_dir = "meshes/"

k = 16

final_results = Matrix{Any}(undef,length(files),6)

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
    prefix = "ex2_" * replace(file_name, ".mat" => "")

    # Run the benchmark function
    println("Processing file: $file_path with prefix: $prefix")
    results = benchmark(A, coords, k, prefix)

    final_results[i,1] = file_name
    final_results[i,2] = results[1]
    final_results[i,3] = results[2]
    final_results[i,4] = results[3]
    final_results[i,5] = results[4]
    final_results[i,6] = results[5]

    # Close the .mat file
    MAT.close(file)
end

# Pretty print the final results
header = [
    "Mesh Name",
    "Recursive($k):part_coordinate",
    "Recursive($k):part_inertial",
    "Recursive($k):part_spectral",
    "part_metis_rec($k)",
    "part_metis_kway($k)"
]
pretty_table(final_results, header=header)