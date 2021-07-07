using Peacock
using PyPlot, LinearAlgebra, LaTeXStrings
r = 0.2
epf(x, y) = x^2 + y^2 < r^2 ? 9 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.01; d2 = 0.01
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, 11, 11)
ky = 0.0

omega0 = 0.01:0.01:1
ii = 1
RT = zeros(length(omega0))
TT = zeros(length(omega0))
for omega in omega0
        kT, dT, VT = solver_EPEW(solver, ky, omega, TM)
        RT[ii], TT[ii] = abcd(kT, dT, VT, omega, ky)
        ii = ii + 1;
end
plot(omega0,RT)