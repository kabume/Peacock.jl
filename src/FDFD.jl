module FDFD

using Peacock

include("solver_fdm.jl")
export Solver, solve, FDFDBasis

end