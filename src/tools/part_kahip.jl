"""
   run_kahip(A::SparseMatrixCSC, np::Int, preconfiguration::String, dbg::Bool=false) -> Tuple{String, Vector{Int}, String}

Run KaHIP patitioning on a graph matrix file.

# Arguments
- `A`: Adjacency matrix of the graph (can be sparse).
- `np`: Number of partitions.
- `preconfiguration`: KaHIP preconfiguration mode.
- `dbg`: Return KaHIP output and file path if true.

# Returns
A tuple containing:
- `stdout::String`: Command output.
- `partitions::Vector{Int}`: Partitioning results.
- `partition_file::String`: Output file path.

# Example
```julia-repl
julia> run_kahip(A, "./graphs/file.graph", 4, "fast")

"""
function run_kahip(A::SparseMatrixCSC, np::Int, preconfiguration::String="fast", dbg::Bool=false)

    mat_path = _write_partition_input(A)

    kahip_exec = _get_executable("KAHIP_PATH")

	cmd = `$kahip_exec $mat_path --k $np --preconfiguration=$preconfiguration`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("KaHIP execution failed with error:\n" * result.code)
    end

    # TODO: Has to be coherent for all external software
    partition_file = "tmppartition" * string(np)
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