using Peacock
using PyPlot
#using Peacock.FDFD

Polarisation = TM
N_eig = 4
eps0 = 1; mu0 = 1
# build geometry
#mu1 = [1 0 0; 0 1 0; 0 0 1]
eps1 = [15 0 0; 0 15 0; 0 0 15]
mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 1]
#eps1 = [14 -12.4im 0; 12.4im 14 0; 0 0 15]
if Polarisation == TE #TODO add to another function
    eps1 = inv(eps1)
else
    mu1 = inv(mu1)
end

r0 = 0.11
function epf(x, y)
    if x^2 + y^2 < r0^2
        return MaterialTensor(eps1)
    else
        return MaterialTensor(eps0)
    end
end

muf(x, y) = MaterialTensor(mu0)

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector

d1 = 0.005 - 1e-6  # resolution along first lattice vector
d2 = 0.005 - 1e-6  # resolution along second lattice vector
#d1 = d2 = 0.2
geometry = Geometry(epf, muf, a1, a2, d1, d2)

#creat ks by Peacock.sample_path
ks = [[[0, 0]]; [[pi, 0]]; [[pi, pi]]; [[2*pi, 2*pi]]]
ks, _ = Peacock.sample_path(ks, dk=2*pi/19)

dx = dy = 1/size(geometry.ep.epsxx)[1]*2
solver_FDFD = Peacock.FDFD.Solver(geometry, [dx dy])
#Peacock.FDFD.Solver(geometry::Geometry)

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

Nx = Int(size(epszz)[1]/2)
Ny = Int(size(epszz)[2]/2)
field_2d = zeros(Nx, Ny)*im

for m = 1:Nx
    for n = 1:Ny
        index = Int((n - 1)*Nx + m)
        field_2d[m, n] = modes[2].data[index]
    end
end