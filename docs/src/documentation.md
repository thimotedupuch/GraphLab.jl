# Documentation
## Graph Construction
```@docs
GraphLab.build_adjacency
```
## Native Partitioning Methods
```@docs
GraphLab.part_coordinate
GraphLab.part_inertial
GraphLab.part_spectral
GraphLab.part_metis
```
## Recursive Bisection
```@docs
GraphLab.recursive_bisection
```
## Visualization
```@docs
GraphLab.draw_graph
```

## Utility Function
```@docs
GraphLab.count_edge_cut
GraphLab.compute_partition_balance
```

## Interfacing with External Partitioning Tools
To integrate with graph partitioning tools such as KaHIP, Graclus, and MDC, we provide utility functions for exporting graph structures, executing external commands, and managing dependencies. These functions streamline the workflow, ensuring compatibility with external software while maintaining interoperability with `GraphLab.jl`.

### Writing a Graph File for Partitioning Tools
This function generates a graph file compatible with tools like KaHIP, Graclus, and MDC from a symmetric adjacency matrix.
**Steps:**
    1. Extract the number of nodes and edges from the adjacency matrix.
    2. Create a temporary file or save it to ./adjacency_lists depending on the temp flag.
    3. Write the graph format:
        - The first line contains the number of nodes and edges.
        - Each subsequent line lists the neighbors of a node.

**Usage:**
    - If temp=true, the file is temporary.
    - If temp=false, the file is saved in ./adjacency_lists with a specified prefix.

```julia
function _write_partition_input(A::SparseMatrixCSC; prefix::String="graph", temp::Bool=false)
    n_nodes = size(A, 1)
    I, J, _ = findnz(A)
    n_edges = div(length(I), 2)

    if temp
        mktemp() do file, _
            write(file, "$n_nodes $n_edges\n")

            for i in 1:n_nodes  
                neighbors = J[findall(==(i), I)]
                write(file, join(neighbors, " ") * "\n")
            end

            return file
        end
    else
        dir = "./adjacency_lists"
        mkpath(dir)
        file_path = joinpath(dir, "$prefix.graph")

        open(file_path, "w") do file
            write(file, "$n_nodes $n_edges\n")
            for i in 1:n_nodes
                neighbors = J[findall(==(i), I)]
                write(file, join(neighbors, " ") * "\n")
            end
        end

        return file_path
    end
end
```

### Running an External Command and Capturing Output
This function runs a shell command and captures both stdout and stderr.
**Steps:**

    1. Execute the command asynchronously.
    2. Capture stdout and stderr into separate tasks.
    3. Wait for the process to finish and return the outputs.


```julia
function _run_cmd(cmd::Cmd)
	out = Pipe()
	err = Pipe()

	process = run(pipeline(cmd, stdout=out, stderr=err), wait=false)

	close(out.in)
	close(err.in)

	stdout_task = @async String(read(out))
	stderr_task = @async String(read(err))

	wait(process)

	return (
		stdout = fetch(stdout_task),
		stderr = fetch(stderr_task),
		code = process.exitcode
	)
end
```

### Retrieving and Validating an Executable Path
This function retrieves an external executable's path from an environment variable and ensures it exists.
**Steps:**

    1. Check if the environment variable is set.
    2. Retrieve the executable path.
    3. Verify that the file exists; otherwise, throw an error.

```julia
function _get_executable(exec_name::String)
    if !haskey(ENV, exec_name)
        error("❌ Environment variable $exec_name is not set. Please define it with the path to the executable.")
    end

    exec_path = ENV[exec_name]

    # Check if the file exists
    if !isfile(exec_path)
        error("❌ Executable '$exec_path' not found. Please check the $exec_name environment variable.")
    end

    return exec_path
end
```
---
With these core utility functions in place, we can now implement functions tailored to specific partitioning software.

### Graclus
Graclus is a spectral-based graph clustering and partitioning software.  

- **Installation**: You need to manually **download and compile Graclus** from its [official source](https://www.cs.utexas.edu/~dml/Software/graclus.html).  
- **Environment Variable**: Set `GRACLUS_PATH` to point to the compiled Graclus executable.  

```bash
export GRACLUS_PATH="/path/to/graclus"
```
```julia
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
```

#### Citation
* Weighted Graph Cuts without Eigenvectors: A Multilevel Approach, I. Dhillon, Y. Guan, and B, Kulis, IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), vol. 29:11, pages 1944-1957, November 2007.
* A Fast Kernel-based Multilevel Algorithm for Graph Clustering, I. Dhillon, Y. Guan, and B, Kulis, Proceedings of The 11th ACM SIGKDD, Chicago, IL, August 21-24, 2005.
* Kernel k-means, Spectral Clustering and Normalized Cuts, I. Dhillon, Y. Guan, and B. Kulis, Proceedings of The 10th ACM SIGKDD, Seattle, WA, August 22-25, 2004.
---

### KaHIP `(run_kahip)`
KaHIP (Karlsruhe High-Quality Partitioning) is an efficient graph partitioning tool that supports various heuristics and optimizations.

- **Installation**: You need to **download and compile KaHIP** from its [official source](https://kahip.github.io/).
- **Environment Variable**: Set `KAHIP_PATH` to the KaHIP binary.

```bash
export KAHIP_PATH="/path/to/kahip"
```
```julia
function run_kahip(A::SparseMatrixCSC, np::Int; preconfiguration::String="fast", dbg::Bool=false)

    mat_path = _write_partition_input(A)

    kahip_exec = _get_executable("KAHIP_PATH")

	cmd = `$kahip_exec $mat_path --k $np --preconfiguration=$preconfiguration`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("KaHIP execution failed with error:\n" * result.code)
    end

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
```
#### Citation
* Peter Sanders and Christian Schulz. Engineering Multilevel Graph Partitioning Algorithms. In Proceedings of the 19th European Symposium on Algorithms (ESA'11), volume 6942 of LNCS, pages 469--480. Springer, 2011.
---
### Multi-Level Diffusion Clustering (MDC) `(run_mdc)`
MDC is an experimental diffusion-based graph partitioning method.

- **Installation**: You must **download, compile, and install MDC** following the instructions from [official source](https://kahip.github.io/).
- **Environment Variable**: Set `MDC_PATH` to point to the MDC executable. It also requires setting `LD_LIBRARY_PATH` for proper execution.
```bash
export MDC_PATH="/path/to/mdc"
export LD_LIBRARY_PATH="/path/to/stag/stag-1.3.0/build_dir/stag_lib/"
```
```julia
function run_mdc(A::SparseMatrixCSC,
                np::Int;
                cut_type::String="ncut",
                local_search::Int=0,
                beta_testing::Bool=true,
                spectral_method::Bool=true,
                beta_min::Int=0,
                beta_max::Int=2,
                dbg::Bool=false)

    mat_path = _write_partition_input(A)

    mdc_exec = _get_executable("MDC_PATH")

	cmd = `$mdc_exec $mat_path $np -b $beta_testing -s $spectral_method -O $cut_type -y $beta_min -z $beta_max -L $local_search`;

    result = _run_cmd(cmd)
	
    if result.code != 0
        error("❌ MDC execution failed with error:\n" * result.code)
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
```
#### Citation
* Lechekhab M., Pasadakis D., Schenk O. (forthcoming) Multilevel Diffusion Based Spectral Graph Clustering. 28th Annual IEEE High Performance Extreme Computing Virtual Conference. IEEE. 28th Annual IEEE High Performance Extreme Computing Virtual Conference. Virtual Conference. September 23-27, 2024

