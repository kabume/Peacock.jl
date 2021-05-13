using Arpack,SparseArrays
using Peacock
using PyPlot
get_k = Peacock.get_k
include("simple-2d-fdfd-master/fdfd/diff_yee2.jl")
"""
    EM, f = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc=[0,0])
solve_fdm is an FDM EM solver on a Yee grid. It returns WF, the Hz or Ez components of the wavefunctions and f their frequency in inverse unit of length.
- `dx` and `dy` are spatial resolution of the single grid (not the double grid).
- `f0` is the initial guess for frequency.
- `eps2` is the dielectic on a double grid, the real simulated grid has two times less resolution.
- `mu2` is the permeability on a double grid, the real simulated grid has two times less resolution.
- `N_eig` is the number of eigenvalues searched.
- `TE=true` if the mode is TE and false if the mode is TM.
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


TE = true #通过TE/TM来判断取eps和mu原始值还是inv
eps0 = 1; mu0 = 1
mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 15]
eps1 = [14 -12.4im 0; 12.4im 14 0; 0 0 15]
if TE == true
    eps1 = inv(eps1)
else
    mu1 = inv(mu1)
end
struct Mu
    muxx::Matrix{ComplexF64}
    muyy::Matrix{ComplexF64}
    muzz::Matrix{ComplexF64}
    muxy::Matrix{ComplexF64}
    muyx::Matrix{ComplexF64}
end
struct Eps
    epsxx::Matrix{ComplexF64}
    epsyy::Matrix{ComplexF64}
    epszz::Matrix{ComplexF64}
    epsxy::Matrix{ComplexF64}
    epsyx::Matrix{ComplexF64}
end

r0 = 0.11
function epfxx(x, y)
    if x^2 + y^2 < r0^2
        return eps1[1, 1]
    else
        return 1
    end
end

function epfzz(x, y)
    if x^2 + y^2 < r0^2
        return eps1[3, 3]
    else
        return 1
    end
end

function epfxy(x, y)
    if x^2 + y^2 < r0^2
        return eps1[1, 2]
    else
        return 0
    end
end

function epfyx(x, y)
    if x^2 + y^2 < r0^2
        return eps1[2, 1]
    else
        return 0
    end
end

function mufxx(x,y)
    if x^2 + y^2 < r0^2
        return mu1[1, 1]
    else
        return 1
    end
end

function mufzz(x,y)
    if x^2 + y^2 < r0^2
        return mu1[3, 3]
    else
        return 1
    end
end

function mufxy(x,y)
    if x^2 + y^2 < r0^2
        return mu1[1, 2]
    else
        return 0
    end
end

function mufyx(x,y)
    if x^2 + y^2 < r0^2
        return mu1[2, 1]
    else
        return 0
    end
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector

d1 = 0.005 - 1e-6  # resolution along first lattice vector
d2 = 0.005 - 1e-6  # resolution along second lattice vector
#d1 = d2 = 0.2
geometry = Geometry(epfxx, mufxx, a1, a2, d1, d2); muyy = muxx = geometry.mu; epsxx = epsyy = geometry.ep
geometry = Geometry(epfzz, mufzz, a1, a2, d1, d2); muzz = geometry.mu; epszz = geometry.ep 
geometry = Geometry(epfxy, mufxy, a1, a2, d1, d2); muxy = geometry.mu; epsxy = geometry.ep 
geometry = Geometry(epfyx, mufyx, a1, a2, d1, d2); muyx = geometry.mu; epsyx = geometry.ep

size(epszz)
mu2 = Mu(muxx, muyy, muzz, muxy, muyx)
eps2 = Eps(epsxx, epsyy, epszz, epsxy, epsyx)




a = Px = Py = 1
N_sam = 10
f0 = 0.1
N_eig = 3
kx1 = range(0.01, pi / Px,length = N_sam)
ky1 = zeros(N_sam)
kx2 = ones(N_sam) * pi / Px
ky2 = range(0, pi / Py, length = N_sam)
radius = sqrt((pi / Px)^2 + (pi / Py)^2)
radius_sam = range(radius, 0, length = round(Int64, N_sam * sqrt(2)))
kx3 = radius_sam * cos(pi / 4)
ky3 = radius_sam * sin(pi / 4)
kx = vcat(kx1[1:end], kx2[2:end], kx3[2:end - 1], kx1[1])
ky = vcat(ky1[1:end], ky2[2:end], ky3[2:end - 1], ky1[1])
ks = [[[0, 0]]; [[pi, 0]]; [[pi, pi]]; [[2*pi, 2*pi]]]
ks, _ = Peacock.sample_path(ks, dk=2*pi/19)


f = zeros(ComplexF64,N_eig,length(ks))
dx = dy = 1/size(eps2.epszz)[1]*2

for j = 1:length(ks)
    global f
    #kinc = [kx[j] ky[j]]
    kinc = ks[j]
    WF, f[:,j] = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc) 
end

figure()
for m = 1:N_eig
    plot(real.(f[m,:]),"-bo")
end
title("Band structure")
display(gcf())