using Peacock, PyPlot
polar = TE
lambda1 = 0.3
lambda2 = 0.06
lambda3 = 0.0
lambda4 = 0.006
eps0 = 1.5 
a1 = [1, 0]
a2 = [0, 1]
d1, d2 = 0.005, 0.005
kx = 2pi; ky = 2pi

function epf(x, y)
    Xi = lambda4 * eps0*((cos(kx*x) + cos(ky*y)) + (cos(2kx*x) + cos(2ky*y)))
    eps1 = eps0 * (1 + lambda1*(cos(kx*x) + cos(ky*y)) + lambda2*(cos(2kx*x) + cos(2ky*y)) + lambda3*(sin(kx*x)+sin(ky*y)))
    return [eps1 Xi*im 0;-Xi*im eps1 0;0 0 1]
end
muf(x, y) = 1

geometry = Geometry(epf, muf, a1, a2, d1, d2, polar)
solver = Peacock.FDFD.Solver(geometry, 2*[d1,d2])
G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [M,G,X,M]
figure()
plot_band_diagram(solver, ks, TE, dk=0.2, bands = 1:5, frequency_scale=1/2pi)

Chern_number(epf, muf, a1, a2, d1, d2, polar, Nkx = 8, Nky = 8, bands=1:5)
