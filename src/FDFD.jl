module FDFD

using LaTeXStrings
using SparseArrays
using Arpack
using Peacock

include("solve_fdm.jl")
export Solver, solve

end