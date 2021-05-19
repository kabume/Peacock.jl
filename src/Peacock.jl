module Peacock

using PyPlot
using LinearAlgebra, FFTW
using CUDA

include("utils.jl")
export diff_yee2

include("geometry.jl")
export Geometry, Mu, Eps, MaterialTensor

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate

include("modes.jl")
export Mode, HilbertSpace, get_field, Mode_FDFD, get_field_FDFD, FDFDBasis

include("solver.jl")
export Solver, TE, TM, solve

include("plotting.jl")
export plot

include("FDFD.jl")

include("band_diagrams.jl")
export plot_band_diagram

include("wilson_loops.jl")
export plot_wilson_loop_winding

# Submodule that defines my commonly used crystals
include("Zoo.jl")

end # module
