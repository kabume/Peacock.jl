module FDFD

using Arpack, SparseArrays
using LaTeXStrings

include("geometry.jl")
export Geometry, Mu, Eps, MaterialTensor

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate

include("modes.jl")
export Mode_FDFD, get_field_FDFD

include("utils.jl")

include("solve_fdm.jl")
export Solver, TE, TM, solve

end