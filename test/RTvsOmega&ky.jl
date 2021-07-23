using Peacock
using PyPlot, LinearAlgebra
p = TE
cutoff = 7
#epf(x, y) = x^2 + y^2 < r^2 ? 11.43 : 1
r0 = 0.3; d = 0.32
epf(x, y) = abs(x+d) <=r0/2  && abs(y+d) <= r0/2 || abs(x-d) <=r0/2  && abs(y+d) <= r0/2 || abs(x+d) <=r0/2  && abs(y-d) <= r0/2 || abs(x-d) <=r0/2  && abs(y-d) <= r0/2 ? 1 : 25
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.01*100/121; d2 = 0.01*100/121
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, cutoff, cutoff)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
figure(); plot_band_diagram(solver, ks, p, bands=1:3, frequency_scale=1/2pi)
grid()

#compute transmission and reflection
kym = (-1:0.1:1)*0.5
omega0 = 0.005:0.005:0.3
RT = zeros(ComplexF64, length(kym), length(omega0))
TT = zeros(ComplexF64, length(kym), length(omega0))
jj = 1
for ky in kym
    ii = 1
    for omega in omega0
        kT, dT, VT = solver_EPWE(solver, ky, omega, p)
        RT[jj, ii], TT[jj, ii], _, _ = abcd(kT, dT, VT, omega, ky)
        ii = ii + 1
    end
    jj = jj +1
end

kymm = repeat(kym, 1, length(omega0))
omega00 = repeat(omega0, 1, length(kym))

pcolor(kymm, transpose(omega00),abs.(RT))