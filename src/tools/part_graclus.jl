"""
    run_graclus(A::SparseMatrixCSC, np::Int, cut_type::String="ncut", local_search::Int=0, dbg::Bool=false)

Partition the graph `A` into `np` parts using Graclus.

# Arguments
- `A`: Adjacency matrix (SparseMatrixCSC).
- `np`: Number of partitions.
- `cut_type`: Partitioning objective (`"ncut"` or `"rcut"`, default: `"ncut"`).
- `local_search`: Local search level (default: `0`).
- `dbg`: Enable debug output (default: `false`).

# Returns
- Partitioning result or error if the process fails.
"""
function run_graclus(A::SparseMatrixCSC, np::Int, cut_type::String="ncut", local_search::Int=0, dbg::Bool=false)

    mat_path = _write_partition_input(A)

    graclus_exec = _get_executable("GRACLUS_PATH")

	cmd = `$graclus_exec $mat_path $np -O $cut_type -L $local_search`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("❌ Graclus execution failed with error:\n" * result.code)
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
        println("✅ Partitions extracted from ", partition_file)
    catch e
        error("❌ Failed to read partition file '$partition_file': ", e)
    end

    if dbg
        return result.stdout, partitions, partition_file
    else
        return partitions
    end
end