using Peacock
using PyPlot

Polarisation = Peacock.FDFD.TE
N_eig = 4
# build geometry
eps0 = 1; mu0 = 1
mu1 = [14 0 0; 0 14 0; 0 0 14]
eps1 = [15 0 0; 0 15 0; 0 0 15]
#mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 15]
#eps1 = [14 -12.4im 0; 12.4im 14 0; 0 0 15]
if Polarisation == Peacock.FDFD.TE
    eps1 = inv(eps1)
else
    mu1 = inv(mu1)
end

r0 = 0.15; d = 0.34; w = 0; l = 0.05
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
mu2 = Peacock.FDFD.Mu(muxx, muyy, muzz, muxy, muyx)
eps2 = Peacock.FDFD.Eps(epsxx, epsyy, epszz, epsxy, epsyx)

#creat ks by Peacock.sample_path
ks = [[[0, 0]]; [[pi, 0]]; [[pi, pi]]; [[2*pi, 2*pi]]]
ks, _ = Peacock.sample_path(ks, dk=2*pi/19)

dx = dy = 1/size(eps2.epszz)[1]*2
solver_FDFD = Peacock.FDFD.Solver(eps2, mu2, dx, dy)

function my_solve(k)
    modes = Peacock.FDFD.solve(solver_FDFD,k,Polarisation,bands=1:N_eig)
    return [mode.frequency for mode in modes]
end

f = zeros(ComplexF64,N_eig,length(ks))
for i = 1:length(ks)
    f[:,i] = my_solve(ks[i])[1:N_eig]/2/pi
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