using Graphs, AMD, GraphsMatching, JuMP, Cbc, GraphsOptim

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

function get_maximum_matching(g::Graph)
  match = maximum_weight_matching(g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0, MOI.Silent() => true));
  # match.weight ≈ 1

  return match
end

"""
    compute_vertex_separator(p1border, p2border, gborder, match)

Given:
- `gborder`: a bipartite graph representing the interface (border) between two partitions.
- `p1border`: list of global node indices on the left side of the bipartition (partition A).
- `p2border`: list of global node indices on the right side (partition B).
- `match.mate`: array representing a maximum matching on `gborder`, where `match.mate[i] = j`
  means local node `i` is matched to `j`, and `-1` means unmatched.

This function implements the **constructive proof of König's theorem** (without using flow algorithms) to extract a **minimum vertex cover** in the bipartite border graph, which forms a **vertex separator**.

### Theory:
By König’s theorem, in a bipartite graph, the size of a **maximum matching** equals the size of a **minimum vertex cover**.

This routine uses the standard alternating-path-based construction:
1. Let `U` be the set of unmatched vertices on side `A` (`p1border`)
2. Let `Z` be the set of vertices reachable from `U` by alternating paths
   (paths alternating between non-matching and matching edges)
3. Then, the minimum vertex cover is: K = (A / Z) ∪ (B ∩ Z)
That is, take all unvisited matched vertices on the left side, and all visited vertices on the right side.

### Returns:
- `sep`: the full separator set (union of both sides)
- `p1nonsep`: nodes in `p1border` not in the separator
- `p2nonsep`: nodes in `p2border` not in the separator

This separator has minimal cardinality among all vertex sets that cover all cut edges between the partitions.
"""
function _vtxsep(gborder, p1border, p2border, match)
  offset = length(p1border)
  total = nv(gborder)
  global_border_nodes = vcat(p1border, p2border)

  visited = falses(total)
  queue = Vector{Tuple{Int, Bool}}()  # (node, phase): false = non-matching edge, true = matching edge

  # Start from unmatched p1 nodes (partition A), phase=false (expect non-matching)
  for i in 1:offset
      if match.mate[i] == -1
          push!(queue, (i, false))
          visited[i] = true
      end
  end

  while !isempty(queue)
      (u, phase) = popfirst!(queue)
      for v in neighbors(gborder, u)
          if !visited[v]
              if phase == false
                  if match.mate[u] != v  # non-matching edge
                      visited[v] = true
                      push!(queue, (v, true))
                  end
              else
                  if match.mate[u] == v  # matching edge
                      visited[v] = true
                      push!(queue, (v, false))
                  end
              end
          end
      end
  end

  # Build separator using the constructive proof
  sepp1 = Set{Int}()
  sepp2 = Set{Int}()

  for i in 1:offset
      if match.mate[i] != -1 && !visited[i]
          push!(sepp1, p1border[i])
      end
  end

  for j in 1:(total - offset)
      v = offset + j
      if visited[v]
          push!(sepp2, p2border[j])
      end
  end

  sep = collect(union(sepp1, sepp2))
  return sep, setdiff(p1border, sepp1), setdiff(p2border, sepp2)
  
end


function nested_dissection(A::SparseMatrixCSC, method::Function; coords::Union{Matrix, Nothing}=nothing, minsep::Int=5, verbose::Bool=false)
    n = size(A, 1)
    p = Vector{Int}(undef, n)

    pA = _nested_dissection(A, method; coords = coords, minsep = minsep, verbose = verbose)

    return pA
end


function _nested_dissection(A::SparseMatrixCSC, method::Function; coords::Union{Matrix, Nothing}=nothing, minsep::Int, verbose::Bool)
    n = size(A, 1)
    perm = zeros(Int, n)

    # Find connected components
    g = SimpleGraph(A)
    components = connected_components(g)
    comp_count = size(components)[1]
    offset = 0

   # @info "🔪 Number of connected components: $comp_count"
    numbered = 0
    for (i, component) in enumerate(components)

      nC = length(component)
      A_sub = A[component, component]
      # @info "🧩 Component $i with $nC nodes"
      if nC <= minsep
      # @info "🍃 Leaf reached with $nC nodes"
        pA = symamd(SparseMatrixCSC{Float64, Int64}(A_sub)) # Approximate Minimum Degree ordering
      else
        sub_coords = coords === nothing ? nothing : coords[component, :]
        if sub_coords === nothing
            p_sub = method(A_sub)
        else
            p_sub = method(A_sub, sub_coords)
        end
        # p_sub = method(A_sub, sub_coords)
        p_sub1 = findall(==(1), p_sub)
        p_sub2 = findall(==(2), p_sub)

            # Find border nodes between p1 and p2
            Abord = A_sub[p_sub1, p_sub2]
            rows, cols, _ = findnz(Abord)
            p1border = unique(p_sub1[rows])
            p2border = unique(p_sub2[cols])

            # if isempty(p1border)
            # # Disconnected, sperator is empty
            # @warn "❌ Separator is empty (disconnected subgraph?)"

            # return Int[], p_sub1, p_sub2, p_sub
            # end

            offset = length(p1border)

        # Build the graph of border vertices
        gborder = Graph(offset + length(p2border))

        for i in 1:offset
            for j in 1:length(p2border)
            
            if A_sub[p1border[i], p2border[j]] != 0
                add_edge!(gborder, i, offset + j)
            end
            end
        end

        # Compute maximum matching
        match =  get_maximum_matching(gborder)
        # @show match.mate
        sep,_,_ = _vtxsep(gborder, p1border, p2border, match);
        
        # TODO: Decide if we want to keep this out of the box
        # min vertex cover impl. or if we want ours (above).
        
        # global_border_nodes = vcat(p1border,p2border)
        # local_sep = min_vertex_cover(gborder; integer=false, optimizer=optimizer_with_attributes(Cbc.Optimizer))
        # sep = global_border_nodes[local_sep]

        vtx1 = setdiff(p_sub1, sep)
        vtx2 = setdiff(p_sub2, sep)
        
        if verbose
          display(draw_graph(A_sub, sub_coords, p_sub))
          wait_for_key("Press Enter to continue...")
        end
        # @show(sep)
        # @show(vtx1)
        # @show(vtx2)
        
        # @info "↙️ Recursing on vtx1 (size = $(length(vtx1)))"
        q1 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx1, vtx1], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx1, vtx1], method; coords=sub_coords[vtx1, :], minsep=minsep, verbose=verbose)

        # @info "↘️ Recursing on vtx2 (size = $(length(vtx2)))"
        q2 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx2, vtx2], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx2, vtx2], method; coords=sub_coords[vtx2, :], minsep=minsep, verbose=verbose)
        
        if q1 === nothing
          q1 = []
        else
          # @show q1
          q1 = vec(Int.(q1))
          # q1 = filter(!=(0), q1)
        end
    
        if q2 === nothing
          q2 = []
        else
          # @show q2
          q2 = vec(Int.(q2))
          # q2 = filter(!=(0), q2)
        end

        pA = vcat(vtx1[q1], vtx2[q2], sep)
        # @info "🧩 Assembled perm length = $(length(pA)), expected = $(size(A_sub, 1))"


        if length(pA) != size(A_sub, 1)
            @warn "Wrong permutation length" length_pA = length(pA) expected = size(A_sub, 1)
            wait_for_key("Press Enter to continue...")
            return collect(1:size(A_sub, 1))  # fallback: identity permutation
        end


        # @show(pA)
        #return pA
        end
        perm[numbered + 1:numbered + nC] = component[pA]
        numbered = numbered + nC

        # @info "✅ Final permutation size: $(length(perm))"

    end
    return perm
end