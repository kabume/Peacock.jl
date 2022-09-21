module Peacock

using Base: Number
using LinearAlgebra: include
using PyPlot
using LinearAlgebra, FFTW
#using CUDA

include("utils.jl")
export diff_yee2, eye, diag_R, mask

include("geometry.jl")
export Geometry, Mu, Eps, MaterialTensor, TE, TM, Polarisation

include("plane_wave_basis.jl")
export BrillouinZoneCoordinate, FDFDBasis

include("modes.jl")
export Mode, HilbertSpace, get_field, Mode_FDFD, get_field_FDFD, HilbertSpace_FDFD

include("solver.jl")
export Solver, solve

include("plotting.jl")
export plot

include("FDFD.jl")

include("band_diagrams.jl")
export plot_band_diagram

include("wilson_loops.jl")
export plot_wilson_loop_winding

# Submodule that defines my commonly used crystals
include("Zoo.jl")

include("Chern_number.jl")
export Chern_number

include("solver_EPWE.jl")
export solver_EPWE, abcd

end # module
