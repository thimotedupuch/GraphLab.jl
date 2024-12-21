# M.L. for High Performance Computing Lab @USI & @ETHZ - malik.lechekhab@usi.ch 

"""
    rec_bisection(method, levels, A, coords=zeros(0), vn=zeros(0))

Compute recursive partitioning of graph `A` using a specified `method` and
number of `levels`.

If the `method` is `coords`-based, coordinates must be passed.

# Examples
```julia-repl
julia> rec_bisection("spectralPart", 3, A)
 5
 ⋮
 2

julia> rec_bisection("coordinatePart", 3, A, coords)
 1
 ⋮
 8
```
"""

function _recursive_bisection(method::Function, levels::Int, A::AbstractSparseMatrix, coords::Union{Matrix, Nothing}=nothing, minpoints::Int=8,vn::Vector{Int}=Int[])

    n = size(A)[1]

    if isempty(vn)
        vn = collect(1:n)
    end

    if n < minpoints || levels < 1
        return ones(Int, n)
    else
        if coords !== nothing
            p = method(A, coords)
            idx1 = findall(x -> x == 1, p)
            idx2 = findall(x -> x == 2, p)
            coords1 = coords[idx1,:]
            coords2 = coords[idx2,:]
        else
            p = method(A)
            idx1 = findall(x -> x == 1, p)
            idx2 = findall(x -> x == 2, p)
            coords1 = coords2 = nothing
        end
        
        vn1 = vn[idx1]
        vn2 = vn[idx2]

        A1 = A[idx1, idx1]
        A2 = A[idx2, idx2]

        # if !isemtpy(coords)
        if levels > 1
            levels = levels - 1
            p1 = _recursive_bisection(method, levels, A1, coords1,minpoints, vn1)
            p2 = _recursive_bisection(method, levels, A2, coords2,minpoints, vn2)

            return vcat(p1, p2.+maximum(p1))[sortperm(vcat(vn1, vn2))]
            
        end

        return p[sortperm(vn)]
        # end

    end
end

function recursive_bisection(method::Function, k::Int, A::AbstractSparseMatrix, coords::Union{Matrix, Nothing}=nothing, minpoints::Int=8)

    levels = log2(k)
    if !isinteger(levels)
        @warn "log2($k) is not an integer. Rounding up to the closest integer."
    end
    
    levels = ceil(Int, levels)

    # we will have 2^levels number of partitions. 
    return _recursive_bisection(method,levels,A,coords,minpoints)
end
