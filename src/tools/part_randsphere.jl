"""
Stereographically lift points to the unit sphere in (d+1)-dimensional space.

coords: n × d matrix of points
Returns: n × (d+1) matrix on the unit sphere
"""
function stereoup(coords::Matrix{Float64})
    n, d = size(coords)
    lifted = zeros(n, d+1)

    for i in 1:n
        x = coords[i, :]
        r2 = dot(x, x)
        lifted[i, 1:d] .= 2 * x
        lifted[i, d+1] = r2 - 1
        lifted[i, :] ./= r2 + 1
    end

    return lifted
end


"""
Approximate centerpoint by coordinate-wise median of a random sample.

Arguments:
- X::Matrix{Float64}: n × d data matrix
- sample_size::Int: number of random rows to sample

Returns:
- center::Vector{Float64}: estimated centerpoint
"""
function centerpoint(X::Matrix{Float64}, sample_size::Int)
    n, d = size(X)
    idx = rand(1:n, min(sample_size, n))
    sample = X[idx, :]
    return mapslices(median, sample; dims=1)[:]
end


"""
Conformal map that moves `center` to the origin on the unit sphere.

Arguments:
- center::Vector{Float64}: centerpoint on the sphere (length d+1)
- X::Matrix{Float64}: n × (d+1) matrix of lifted points

Returns:
- Y::Matrix{Float64}: transformed points (same size as X)
"""
function conmap(c::Vector{Float64}, xyz::Matrix{Float64})
    d = size(xyz, 2)
  
    # Compute reflection and stretch
    Q, r = reflector(c)
    alpha = sqrt((1+r)/(1-r))
  
    # Reflect
    xyzref = xyz * Q;
  
    # Handle north pole in the stereographic projection
    norths = findall(x -> abs(x - 1) < sqrt(eps()), xyzref[:, d]) 
    if isempty(norths)
      xyzn = xyzref[norths, :]
      xyzref[norths, :] = zeros(length(norths), d)
    end
  
    xyref = stereodown(xyzref)
  
    # Stretch
    xymap = xyref ./ alpha
    xyzmap = stereoup(xymap)
  
    if isempty(norths)
      xyzmap[norths,:] = xyzn
    end
  
    return xyzmap, xymap
end

function reflector(c)
    d = length(c)
    p = collect(d:-1:1)
    Q, r = qr(c[collect(d:-1:1)])
    Q = Q[p, p]
    r = r[1]
    return Q, r
end

function stereodown(xyz)
    n, dim = size(xyz)
    xy = xyz[:, 1:dim-1] ./ (1 .- xyz[:, dim])  # broadcasting does the repetition
  
    return xy
end


"""
Try ntrials random great circle cuts on unit-sphere data.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- X::Matrix{Float64}: n × (d+1) unit-sphere coordinates
- ntrials::Int: number of random directions to try

Returns:
- bestdir::Vector{Float64}: direction vector of best cut
- mincut::Float64: minimum edge cut found
"""
function sepcircle(A::SparseMatrixCSC, X::Matrix{Float64}, ntrials::Int)
    _, d = size(X)
    M = (X' * X) ^ 2
  
    vv = randn(ntrials, d)*M
    
    bestcut = Inf;
    bestdir = vv[1, :]
  
    for i in 1:size(vv, 1)
      v = vv[i, :]
  
      if norm(v) != 0
        currentcut = sepquality(v, A, X)
      end
  
      if bestcut > currentcut
        bestcut = currentcut
        bestdir = v
      end
    end
  
    return bestdir, bestcut 
  
end

function sepquality(v::Vector{Float64}, A::SparseMatrixCSC, xyz::Matrix{Float64})
    a, b = partition(xyz, v)
    cutsize = nnz(Bool.(A[a, b]) .| Bool.(A[b, a]'))
  
    return cutsize
end


"""
Try ntrials random straight hyperplanes in Euclidean space.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- coords::Matrix{Float64}: n × d original coordinates
- ntrials::Int: number of directions to try

Returns:
- bestdir::Vector{Float64}: best direction
- mincut::Float64: corresponding edge cut
"""
function sepline(A::SparseMatrixCSC, xy::Matrix{Float64}, ntrials::Int)
    d = size(xy, 2)
  
    _, S, V = svd(xy);
  
    exponent = 2*(d+1)/(ntrials-1)
    s = diagm(S).^exponent
    W = V * s * V'
    vv = randn(ntrials, d) * W
    rownorms = sqrt.(sum(vv .* vv, dims=2))
    vv = Diagonal(vec(1 ./ rownorms)) * vv
  
    bestcut = Inf
    bestdir = vv[1, :]
    for i in 1:ntrials
      v = vv[i, :]
      cut = sepquality(v, A, xy)
      if cut < bestcut
        bestcut = cut
        bestdir = v
      end
    end
  
    return bestdir, bestcut 
  end


"""
Geometric partitioning using random spheres and lines.

Arguments:
- A::SparseMatrixCSC: adjacency matrix
- coords::Matrix{Float64}: n × d node positions
- ntries::Int (optional): number of directions to try (default: 30)

Returns:
- part1, part2: vectors of node indices
"""
function part_randsphere(A::SparseMatrixCSC, coords::Matrix{Float64}; ntrials::Int=30)
    n, d = size(coords)

  # How to split the tries
  nlines = floor(Int, (ntrials / 2) ^ (d / (d + 1))) # number of lines to try (fallback)
  nouter = ceil(Int, log(ntrials - nlines + 1) / log(20)) # number of centerpoints to try
  ninner = floor(Int, (ntrials - nlines) / nouter) # number of random directions (great circles) to test for each centerpoint.
  nlines = ntrials - nouter * ninner # Rounding

  # Normalize coordinates: center and scale
  xy_shift = mean(coords, dims=1)
  xy = coords .- xy_shift
  scale = maximum(abs.(xy))
  xy ./= scale

  # Project to sphere
  xyz = stereoup(xy)

  # Spherical sperator search


  csample = min(n, (d + 3)^4)
  bestcut = Inf
  bestdir = bestcpt = fill(NaN, d + 1)

  for _ in 1:nouter
    cpt = centerpoint(xyz, csample)
    if 1 - norm(cpt) < sqrt(eps())
      dc = randn(size(cpt))
      cpt = 0.9 * cpt + dc / (20 * norm(dc))
    end
    xyzmap, _ = conmap(cpt, xyz) # THINGS WORKING TILL HERE
    circledir, circlecut = sepcircle(A, xyzmap, ninner)
    if circlecut < bestcut
      bestcut = circlecut
      bestdir = circledir
      bestcpt = cpt
    end
  end

  linedir, linecut = sepline(A, xy, nlines)

  if linecut < bestcut
    bestcut = linecut
    bestdir = linedir
    p1, p2 = partition(xy, linedir)
  else
    bestmap, _  = conmap(bestcpt, xyz)
    p1, p2 = partition(bestmap, bestdir)
  end

  p = ones(Int, n)
  p[p1] .= 1
  p[p2] .= 2

  return p
end

function partition(X::Matrix{Float64}, direction::Vector{Float64})
    proj = X * direction
    threshold = median(proj)
    part1 = findall(proj .<= threshold)
    part2 = findall(proj .> threshold)
    return part1, part2
end
