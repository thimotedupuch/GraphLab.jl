module GraphLabExamples
using Pkg.Artifacts, MAT, SparseArrays
export airfoil

"Return (A::SparseMatrixCSC, coords::Matrix) from the airfoil artifact."
function airfoil()
    dir  = artifact"airfoil1"
    file = joinpath(dir, "airfoil", "airfoil1.mat")
    d      = MAT.matread(file)
    A      = sparse(d["Problem"]["A"])
    coords = Matrix(d["Problem"]["aux"]["coord"])
    return A, coords
end
end