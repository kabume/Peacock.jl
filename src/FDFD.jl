module FDFD

using Peacock

include("solve_fdm.jl")
export Solver, solve, FDFDBasis

end