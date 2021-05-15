module FDFD

using Arpack, SparseArrays
using LaTeXStrings

include("geometry.jl")
export Geometry

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate

include("modes.jl")
export Mode, HilbertSpace, get_field

include("utils.jl")

include("solve_fdm.jl")
export Solver, TE, TM, solve

end