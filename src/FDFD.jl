module FDFD

using Peacock
using Arpack, SparseArrays
using LaTeXStrings
using LinearAlgebra

include("solve_fdm.jl")
export Solver, solve

end