using Arpack, SparseArrays
using Peacock
using PyPlot
include("src/diff_yee2.jl")
include("src/solve_fdm.jl")

# build geometry
TE = true 
eps0 = 1; mu0 = 1
#mu1 = [1 0 0; 0 1 0; 0 0 1]
#eps1 = [8.9 0 0; 0 8.9 0; 0 0 8.9]
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

r0 = 0.15; d = 0.35; w = 0; l = 0.05
function epfxx(x, y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return eps1[1, 1]
    else
        return 1
    end
end

function epfzz(x, y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return eps1[3, 3]
    else
        return 1
    end
end

function epfxy(x, y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return eps1[1, 2]
    else
        return 0
    end
end

function epfyx(x, y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return eps1[2, 1]
    else
        return 0
    end
end

function mufxx(x,y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return mu1[1, 1]
    else
        return 1
    end
end

function mufzz(x,y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return mu1[3, 3]
    else
        return 1
    end
end

function mufxy(x,y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
        return mu1[1, 2]
    else
        return 0
    end
end

function mufyx(x,y)
    if (x+d)^2 + (y+d)^2 <= r0^2 || (x-d)^2 + (y-d)^2 <= r0^2 || (x+d)^2 + (y-d)^2 <= r0^2 || (x-d)^2 + (y+d)^2 <= r0^2
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

# creat ks by hand
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

#creat ks by Peacock.sample_path
ks = [[[0, 0]]; [[pi, 0]]; [[pi, pi]]; [[2*pi, 2*pi]]]
ks, _ = Peacock.sample_path(ks, dk=2*pi/19)


#f = zeros(ComplexF64,N_eig,length(ks))
f = zeros(ComplexF64,N_eig,length(kx))
dx = dy = 1/size(eps2.epszz)[1]*2

for j = 1:length(ks)
    global f
    #kinc = [kx[j] ky[j]]
    kinc = ks[j]
    WF, f[:,j] = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc) 
end

f = zeros(ComplexF64,N_eig,length(kx))
dx = dy = 1/size(eps2.epszz)[1]*2

for j = 1:length(kx)
    global f
    kinc = [kx[j] ky[j]]
    WF, f[:,j] = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc) 
end

figure()
for m = 1:N_eig
    plot(real.(f[m,:]),"-bo")
end
title("Band structure")
display(gcf())

Nkx = Nky = 4
shifted=0.0
Pkx = (2*pi + shifted)/Px
Pky = (2*pi + shifted)/Py
deltakx = Pkx/Nkx
deltaky = Pky/Nky
NK = Nkx*Nky
for m = 1:Nkx
    for n = 1:NkY
        index = (n-1) * Nkx + m
        kx[index, :] = [(m-1)*deltakx (m)*deltakx (m)*deltakx (m-1)*deltakx]
        ky[index, :] = [(n-1)*deltaky (n-1)*deltaky (n)*deltaky (n)*deltaky]
        WF, f[:,j] = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc) 


    end
end