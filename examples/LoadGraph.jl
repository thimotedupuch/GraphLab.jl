module LoadGraph
using Pkg.Artifacts, MAT, SparseArrays
export airfoil
function airfoil()
    dir = artifact"airfoil1"
    file = joinpath(dir, "airfoil", "airfoil1.mat")
    d = MAT.matread(file)
    return sparse(d["Problem"]["A"]), Matrix(d["Problem"]["aux"]["coord"])
end
end