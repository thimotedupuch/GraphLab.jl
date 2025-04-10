using Graphs, AMD, GraphsMatching, JuMP, Cbc

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

function get_maximum_matching(g::Graph)
  match = maximum_weight_matching(g, optimizer_with_attributes(Cbc.Optimizer, "LogLevel" => 0, MOI.Silent() => true));
  # match.weight ≈ 1

  return match
end


function _vtxsep(gborder, p1border, p2border, match)
    # König’s theorem says: in bipartite graphs, the size of a maximum matching = size of a minimum vertex cover.
    # The minimum vertex cover can be extracted from the alternating tree.
    offset = length(p1border)
    total = nv(gborder)
  
    # 1. Start from unmatched node
    visited = falses(total);
    queue = Int[]
    for i in 1:offset
      if match.mate[i] == -1
        push!(queue, i)
        visited[i] = true
      end
    end
  
    # 2. BFS alternating tree
    while !isempty(queue)
      u = popfirst!(queue)
  
      for v in neighbors(gborder, u)
        if !visited[v]
          if match.mate[u] != v
            visited[v] = true
            push!(queue, v)
          end
        end
      end
    end
  
    # 3. Build separator
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
  
    sep = union(sepp1, sepp2) |> collect
  
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
    perm = zeros(1, n)

    # Find connected components
    g = SimpleGraph(A)
    components = connected_components(g)
    comp_count = size(components)[1]
    offset = 0

   @info "🔪 Number of connected components: $comp_count"
    numbered = 0
    for (i, component) in enumerate(components)

        nC = length(component)
        A_sub = A[component, component]
        @info "🧩 Component $i with $nC nodes"
        if nC <= minsep
        @info "🍃 Leaf reached with $nC nodes"
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

            if isempty(p1border)
            # Disconnected, sperator is empty
            @warn "❌ Separator is empty (disconnected subgraph?)"
            return Int[], p_sub1, p_sub2, p_sub
            end

            offset = length(p1border)

        # Build the graph of border vertices
        gborder = Graph(offset + length(p2border))

        for i in 1:offset
            for j in 1:length(p2border)
            if A[p1border[i], p2border[j]] != 0
                add_edge!(gborder, i, offset + j)
            end
            end
        end

        # Compute maximum matching
        match =  get_maximum_matching(gborder)

        sep,_,_ = _vtxsep(gborder, p1border, p2border, match);

        vtx1 = setdiff(p_sub1, sep)
        vtx2 = setdiff(p_sub2, sep)
        
        if verbose
          display(draw_graph(A_sub, sub_coords, p_sub))
          wait_for_key("Press Enter to continue...")
        end
        # @show(sep)
        # @show(vtx1)
        # @show(vtx2)
        
      @info "↙️ Recursing on vtx1 (size = $(length(vtx1)))"
      q1 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx1, vtx1], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx1, vtx1], method; coords=sub_coords[vtx1, :], minsep=minsep, verbose=verbose)

      @info "↘️ Recursing on vtx2 (size = $(length(vtx2)))"
      q2 = sub_coords === nothing ? 
                _nested_dissection(A_sub[vtx2, vtx2], method; minsep=minsep, verbose=verbose) :
                _nested_dissection(A_sub[vtx2, vtx2], method; coords=sub_coords[vtx2, :], minsep=minsep, verbose=verbose)
        
      if q1 === nothing
        q1 = []
      else
        q1 = vec(Int.(q1))
        q1 = filter(!=(0), q1)
      end
    
      if q2 === nothing
        q2 = []
      else
        q2 = vec(Int.(q2))
        q2 = filter(!=(0), q2)
      end

        pA = vcat(vtx1[q1], vtx2[q2], sep)
        @info "🧩 Assembled perm length = $(length(pA)), expected = $(size(A_sub, 1))"


        if length(pA) != size(A_sub, 1)
            @warn "Wrong permutation length" length_pA = length(pA) expected = size(A_sub, 1)
            return collect(1:size(A_sub, 1))  # fallback: identity permutation
        end



        # return pA
        end
        perm[numbered + 1:numbered + nC] = component[pA]
        numbered = numbered + nC

        @info "✅ Final permutation size: $(length(perm))"

        return vec(Int.(perm))
    end
end