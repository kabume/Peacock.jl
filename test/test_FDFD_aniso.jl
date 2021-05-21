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
#a2 = 0.5*[sqrt(3), 1]

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

modes = Peacock.FDFD.solve(solver_TM, M, TM, bands=1:4)

angles=[60n for n in 0:5]; ep_bg=1; ep_cyl=11.7

R = 1/2.9; d = 2R/3

epf(x,y) = Peacock.Zoo.wu_ep(x, y, angles, R, d, d, ep_bg, ep_cyl)

muf(x,y)=1

d1 = d2 = 1/500

a1 = [1, 0]
a2 = [0.5, sqrt(3)/2]
geometry = Geometry(epf, muf, a1, a2, 1/500, 1/500, TM)

geometry = Geometry(epf, muf, -60, +60, 1/500, 1/500, TM)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
M = BrillouinZoneCoordinate(  0, 1/sqrt(3), "M")
K = BrillouinZoneCoordinate(1/3, 1/3, "K")

ks = [K,G,M]

solver = Peacock.FDFD.Solver(geometry, [2*d1, 2*d2])

plot_band_diagram(solver, ks, TM, color="red",bands=1:6, dk=0.6, frequency_scale=1/2pi)