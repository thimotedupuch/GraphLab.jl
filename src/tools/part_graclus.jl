"""
    run_graclus(A::SparseMatrixCSC, np::Int, cut_type::String="ncut", local_search::Int=0, dbg::Bool=false)

Partition the graph `A` into `np` parts using Graclus.

# Arguments
- `A`: Adjacency matrix (SparseMatrixCSC).
- `np`: Number of partitions.
- `cut_type`: Partitioning objective (`"ncut"` or `"rassoc"`, default: `"ncut"`).
- `local_search`: Local search level (default: `0`).
- `dbg`: Enable debug output (default: `false`).

# Returns
- `partitions::Vector`: Partitioning results if `dbg = false`.
Or a tuple if `dbg = true`:
- `stdout::String`: Command output.
- `partitions::Vector`: Partitioning results.
- `partition_file::String`: Output file path.
"""
function run_graclus(A::SparseMatrixCSC, np::Int; cut_type::String="ncut", local_search::Int=0, dbg::Bool=false)

    mat_path = _write_partition_input(A)

    graclus_exec = _get_executable("GRACLUS_PATH")

	cmd = `$graclus_exec $mat_path $np -O $cut_type -L $local_search`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("❌ Graclus execution failed with error:\n" * result.code)
    end

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
        println("🖥️ Execution Result:\n", result.stdout)
        println("📂 Partition File: ", partition_file)
        println("🔢 Partitions: ", partitions)
        return partitions
    else
        # Remove the partition file after reading
        rm(partition_file, force=true)
        return partitions
    end
end