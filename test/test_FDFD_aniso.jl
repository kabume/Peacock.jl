using Peacock
using PyPlot
#using Peacock.FDFD

Polarisation = TM
N_eig = 4
eps0 = 1; mu0 = 1
# build geometry
#mu1 = 1
eps1 = 15
mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 15]
#eps1 = [14 -12.4im 0; 12.4im 14 0; 0 0 15]

r0 = 0.11
function epf(x, y)
    if Polarisation == TE
        if x^2 + y^2 < r0^2
            return inv(MaterialTensor(eps1))
        else
            return inv(MaterialTensor(eps0))
        end
    else
        if x^2 + y^2 < r0^2
            return MaterialTensor(eps1)
        else
            return MaterialTensor(eps0)
        end
    end
end

#muf(x, y) = MaterialTensor(mu0)

function muf(x, y)
    if Polarisation == TM
        if x^2 + y^2 < r0^2
            return inv(MaterialTensor(mu1))
        else
            return inv(MaterialTensor(mu0))
        end
    else
        if x^2 + y^2 < r0^2
            return MaterialTensor(mu1)
        else
            return MaterialTensor(mu0)
        end
    end
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector

d1 = 0.005  # resolution along first lattice vector
d2 = 0.005  # resolution along second lattice vector
geometry = Geometry(epf, muf, a1, a2, d1, d2)

#creat ks by Peacock.sample_path
ks = [[[0, 0]]; [[pi, 0]]; [[pi, pi]]; [[2*pi, 2*pi]]]
ks, _ = Peacock.sample_path(ks, dk=2*pi/19)

solver_FDFD = Peacock.FDFD.Solver(geometry, [2*d1, 2*d2])

function my_solve(k)
    modes = Peacock.FDFD.solve(solver_FDFD,k,Polarisation,bands=1:N_eig)
    return [mode.frequency for mode in modes]
end

f = zeros(ComplexF64, N_eig, length(ks))
for i = 1:length(ks)
    f[:,i] = my_solve(ks[i])[1:N_eig]/2/pi
end

figure()
for m = 1:N_eig
    plot(real.(f[m,:]),"-bo")
end
title("Band structure")
display(gcf())

Nx = Int(size(geometry.ep.epsxx)[1]/2)
Ny = Int(size(geometry.ep.epsxx)[2]/2)
field_2d = zeros(Nx, Ny)*im
k = [0, 0]
modes = Peacock.FDFD.solve(solver_FDFD,k,Polarisation,bands=1:N_eig)

for m = 1:Nx
    for n = 1:Ny
        index = Int((n - 1)*Nx + m)
        field_2d[m, n] = modes[2].data[index]
    end
end