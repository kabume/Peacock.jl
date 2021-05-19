using Peacock
using PyPlot
#using Peacock.FDFD

# Permittivity
function epf(x,y)
    # equation of a circle with radius 0.2a
    if x^2+y^2 <= 0.2^2
        # dielectric inside the circle
        return 8.9
    else
        # air outside the circle
        return 1
    end
end

# Permeability is unity everywhere
function muf(x,y)
    return 1
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector

d1 = 0.005  # resolution along first lattice vector
d2 = 0.005  # resolution along second lattice vector

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]

figure(figsize=(4,3))
geometry_FDFD = Geometry(epf, muf, a1, a2, d1, d2, TE)
solver_FDFD = Peacock.FDFD.Solver(geometry_FDFD, [2*d1, 2*d2])
plot_band_diagram(solver_FDFD, ks, TE, color="red",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)

geometry_FDFD = Geometry(epf, muf, a1, a2, d1, d2, TM)
solver_FDFD = Peacock.FDFD.Solver(geometry_FDFD, [2*d1, 2*d2])
plot_band_diagram(solver_FDFD, ks, TM, color="blue",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
ylim(0,0.8)

geometry_PWE = Geometry(epf, muf, a1, a2, d1, d2)
solver_PWE =  Solver(geometry_PWE, 7)
figure(figsize=(4,3))
plot_band_diagram(solver_PWE, ks, TE, color="red",
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
plot_band_diagram(solver_PWE, ks, TM, color="blue";
            bands=1:4, dk=0.1, frequency_scale=1/2pi)
ylim(0,0.8)

modes = Peacock.FDFD.solve(solver_FDFD, X, TM)