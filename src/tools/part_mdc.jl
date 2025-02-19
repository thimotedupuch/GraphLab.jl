"""
    run_mdc(A::SparseMatrixCSC, np::Int; cut_type::String="ncut", local_search::Int=0, 
            beta_testing::Bool=true, spectral_method::Bool=true, beta_min::Int=0, 
            beta_max::Int=2, dbg::Bool=false)

Partition the graph `A` into `np` parts using the MDC algorithm.

# Arguments
- `A`: Adjacency matrix (SparseMatrixCSC).
- `np`: Number of partitions.
- `cut_type`: Partitioning objective (`"ncut"` or `"rcut"`, default: `"ncut"`).
- `local_search`: Local refinement level (default: `0`).
- `beta_testing`: Enable beta parameter tuning (default: `true`).
- `spectral_method`: Use spectral partitioning (default: `true`).
- `beta_min`: Minimum beta value (default: `0`).
- `beta_max`: Maximum beta value (default: `2`).
- `dbg`: Enable debug output (default: `false`).

# Returns
- Partitioning result or error if the process fails.
"""
function run_mdc(A::SparseMatrixCSC,
                np::Int,
                cut_type::String="ncut",
                local_search::Int=0,
                beta_testing::Bool=true,
                spectral_method::Bool=true,
                beta_min::Int=0,
                beta_max::Int=2,
                dbg::Bool=false)

    mat_path = _write_partition_input(A)

    # TODO: Needs to be set by user
    # Set the library path to locate the STAG library
    ENV["LD_LIBRARY_PATH"] = "../Graphs/graclus1.2/stag/stag-1.3.0/build_dir/stag_lib/"

    # TODO: Needs to be set by user
	cmd = `../Graphs/graclus1.2/graclus $mat_path $np -b $beta_testing -s $spectral_method -O $cut_type -y $beta_min -z $beta_max -L $local_search`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("MDC execution failed with error:\n" * result.code)
    end

    # TODO: Has to be coherent for all external software
    partition_file = "graph.graph"*".part."*string(np)
    partitions = Int[]

    try
        open(partition_file, "r") do f
            for line in eachline(f)
                push!(partitions, parse(Int, line))
            end
        end
        println("Partitions extracted from ", partition_file)
    catch e
        error("Failed to read partition file '$partition_file': ", e)
    end

    if dbg
        return result.stdout, partitions, partition_file
    else
        return partitions
    end
end