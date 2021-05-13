"""
    EM, f = solve_fdm(dx::Float64, dy::Float64, f0::Float64, eps2::Eps, mu2::Mu, N_eig::Int64, TE::Bool, kinc=[0 0])
solve_fdm is an FDM EM solver on a Yee grid. It returns WF, the Hz or Ez components of the wavefunctions and f their frequency in inverse unit of length.
- `dx` and `dy` are spatial resolution of the single grid (not the double grid).
- `f0` is the initial guess for frequency.
- `eps2` is the dielectic on a double grid, the real simulated grid has two times less resolution.
- `mu2` is the permeability on a double grid, the real simulated grid has two times less resolution.
- `N_eig` is the number of eigenvalues searched.
- `TE = true` if the mode is TE and false if the mode is TM.
- `kinc` is the wavevector associated with Bloch mode boundary conditions, (Periodic boundary conditions with a phase term).
"""

function solve_fdm(dx::Float64, dy::Float64, f0::Float64, eps2::Eps, mu2::Mu, N_eig::Int64, TE::Bool, kinc=[0 0])
    k0 = 2*pi*f0
    BC = [-2 -2] # periodic boundary conditions
    Nx2 = size(eps2.epszz)[1]
    Ny2 = size(eps2.epszz)[2]
    DEX, DEY, DHX, DHY = diff_yee2(([Nx2 Ny2]/2), [dx dy], BC, kinc)
    EM = zeros(ComplexF64, convert(Int64,Nx2 * Ny2 /4), 3, N_eig)
    if TE # TE mode
        epsxx = eps2.epsxx[2:2:Nx2, 1:2:Ny2]
        epsxx_inv = spdiagm(0 => 1 ./ epsxx[:])
        epsxx = spdiagm(0 => epsxx[:])
        epsyy = eps2.epsyy[1:2:Nx2, 2:2:Ny2]
        epsyy_inv = spdiagm(0 => 1 ./ epsyy[:])
        epsyy = spdiagm(0 => epsyy[:])
        epsxy = eps2.epsxy[1:2:Nx2, 2:2:Ny2]
        epsxy_inv = spdiagm(0 => 1 ./ epsxy[:])
        epsxy = spdiagm(0 => epsxy[:])
        epsyx = eps2.epsyx[2:2:Nx2, 1:2:Ny2]
        epsyx_inv = spdiagm(0 => 1 ./ epsyx[:])
        epsyx = spdiagm(0 => epsyx[:])
        muzz = mu2.muzz[2:2:Nx2, 2:2:Ny2]
        muzz_inv = spdiagm(0 => 1 ./ muzz[:])
        muzz = spdiagm(0 => muzz[:])
        H = DHX * (epsyx*DEY - epsyy*DEX) - DHY * (epsxx*DEY-epsxy*DEX)
        ks, WFs = eigs(H, muzz, nev = N_eig, sigma = k0^2) # solver
    else # TM mode
        epszz = eps2.epszz[1:2:Nx2, 1:2:Ny2]
        epszz_inv = spdiagm(0 => 1 ./ epszz[:])
        epszz = spdiagm(0 => epszz[:])
        muxx = mu2.muxx[1:2:Nx2, 2:2:Ny2]
        muxx_inv = spdiagm(0 => 1 ./ muxx[:])
        muxx = spdiagm(0 => muxx[:])
        muyy = mu2.muyy[2:2:Nx2, 1:2:Ny2]
        muyy_inv = spdiagm(0 => 1 ./ muyy[:])
        muyy = spdiagm(0 => muyy[:])
        muxy = mu2.muxy[2:2:Nx2, 1:2:Ny2]
        muxy_inv = spdiagm(0 => 1 ./ muxy[:])
        muxy = spdiagm(0 => muxy[:])
        muyx = mu2.muyx[1:2:Nx2, 2:2:Ny2]
        muyx_inv = spdiagm(0 => 1 ./ muyx[:])
        muyx = spdiagm(0 => muyx[:])
        H = DEX * (muyx*DHY - muyy*DHX) - DEY * (muxx*DHY - muxy*DHX)
        ks, WFs = eigs(H, epszz, nev = N_eig, sigma = k0^2)
    end

    k = sort(real.(sqrt.(ks)))
    order = sortperm(real.(ks))

    eta_0 = 376.73031346177 # free space impedance
    for i=1:length(ks)
        if TE
            EM[:, 1, i] = WFs[:, order[i]]*im/eta_0 #Hz
            EM[:, 2, i] = epsxx_inv*DHY*WFs[:, order[i]]/k0 #Ex it's not right for anisotropic case
            EM[:, 3, i] = -epsyy_inv*DHX*WFs[:, order[i]]/k0 #Ey
        else
            EM[:, 1, i]=WFs[:, order[i]] #Ez
            EM[:, 2, i]=muxx_inv*DHY*WFs[:, order[i]]*im/eta_0/k0 #Hx
            EM[:, 3, i]=-muyy_inv*DHX*WFs[:, order[i]]*im/eta_0/k0 #Hy
        end
    end
    EM=reshape(EM, convert(Int64, Nx2/2), convert(Int64, Ny2/2), 3, N_eig)
    f = k/(2*pi)
    return EM, f
end
