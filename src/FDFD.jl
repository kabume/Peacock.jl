module FDFD

using Arpack, SparseArrays

include("geometry.jl")
export Geometry

include("modes.jl")
export Mode, HilbertSpace, get_field

include("utils.jl")

include("solve_fdm.jl")
export Solver, TE, TM, solve

end