using LaTeXStrings
using SparseArrays
using Arpack


struct Solver
    eps2::Eps
    mu2::Mu
    basis::FDFDBasis
end

function Solver(geometry::Geometry, resolution::Array{<:Real,1}; GPU=false)
    basis = FDFDBasis(geometry.a1, geometry.a2, resolution)
    return Solver(geometry.ep, geometry.mu, basis)
end


function solve(solver::Solver, k::AbstractVecOrMat{<:Real}, polarisation::Peacock.Polarisation; bands=1:5)
    
    eps2, mu2 = solver.eps2, solver.mu2
    dx, dy = solver.basis.resolution[1], solver.basis.resolution[2]
    BC = [-2 -2]
    Nx2 = size(solver.eps2.epszz)[1]
    Ny2 = size(solver.eps2.epszz)[2]
    DEX, DEY, DHX, DHY = diff_yee2(([Nx2 Ny2]/2), [dx dy], BC, solver.basis.a1, solver.basis.a2, k)

    if polarisation == TE
        epsxx = eps2.epsxx[2:2:Nx2, 1:2:Ny2]
        epsxx = spdiagm(0 => epsxx[:])
        epsyy = eps2.epsyy[1:2:Nx2, 2:2:Ny2]
        epsyy = spdiagm(0 => epsyy[:])
        epsxy = eps2.epsxy[1:2:Nx2, 2:2:Ny2]
        epsxy = spdiagm(0 => epsxy[:])
        epsyx = eps2.epsyx[2:2:Nx2, 1:2:Ny2]
        epsyx = spdiagm(0 => epsyx[:])
        muzz = mu2.muzz[2:2:Nx2, 2:2:Ny2]
        weighting = muzz[:]
        muzz = spdiagm(0 => muzz[:])
        LHS = DHX * (epsyx*DEY - epsyy*DEX) - DHY * (epsxx*DEY-epsxy*DEX)
        RHS = muzz
        label = L"H_z"
    elseif polarisation == TM
        epszz = eps2.epszz[1:2:Nx2, 1:2:Ny2]
        weighting = epszz[:]
        epszz = spdiagm(0 => epszz[:])
        muxx = mu2.muxx[1:2:Nx2, 2:2:Ny2]
        muxx = spdiagm(0 => muxx[:])
        muyy = mu2.muyy[2:2:Nx2, 1:2:Ny2]
        muyy = spdiagm(0 => muyy[:])
        muxy = mu2.muxy[2:2:Nx2, 1:2:Ny2]
        muxy = spdiagm(0 => muxy[:])
        muyx = mu2.muyx[1:2:Nx2, 2:2:Ny2]
        muyx = spdiagm(0 => muyx[:])
        LHS = DEX * (muyx*DHY - muyy*DHX) - DEY * (muxx*DHY - muxy*DHX)
        RHS = epszz
        label = L"E_z"
    end

    freqs_squared, modes_data = try
        eigs(LHS, RHS, nev =bands[end], sigma = (2*pi*0.1)^2)
    catch
        eigs(RHS \ LHS, nev = bands[end], sigma = (2*pi*0.1)^2)
    end

    freqs = sqrt.(freqs_squared)
    # Sort by increasing frequency
    idx = sortperm(freqs, by=real)
    freqs = freqs[idx][bands]
    modes_data = modes_data[:,idx][:,bands]

    modes = Mode_FDFD[]
    for i in 1:length(freqs)
        mode = Mode_FDFD(k, freqs[i], modes_data[:,i], weighting, solver.basis, label)
        push!(modes, mode)
    end
    return modes
end

function solve(solver::Solver, x::BrillouinZoneCoordinate, polarisation::Polarisation; bands=:)
    k = Peacock.get_k(x, solver.basis)
    return solve(solver, k, polarisation, bands=bands)
end