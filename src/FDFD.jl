module FDFD

using Arpack, SparseArrays

include("utils.jl")

include("solve_fdm.jl")
export Solver, TE, TM, solve

end