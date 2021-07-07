using Peacock
using PyPlot
#using Peacock.FDFD

function change_cell(kappa, Nkx, Nky, d1, d2, bands)
#kappa = 0im
mu1 = [14 kappa 0;-kappa 14 0;0 0 1]
# Permittivity
function epf(x,y)
    # equation of a circle with radius 0.2a
    s = 1/2
    r = 0.1/sqrt(2)
    if x^2+y^2 <= r^2 || (x+s)^2+(y+s)^2 < r^2 || (x+s)^2+(y-s)^2 < r^2 || (x-s)^2+(y+s)^2 < r^2 || (x-s)^2+(y-s)^2 < r^2
        # dielectric inside the circle
        return 15
    else
        # air outside the circle
        return 1
    end
end
function muf(x,y)
    # equation of a circle with radius 0.2a
    s = 1/2
    r = 0.1/sqrt(2)
    if x^2+y^2 <= r^2 || (x+s)^2+(y+s)^2 < r^2 || (x+s)^2+(y-s)^2 < r^2 || (x-s)^2+(y+s)^2 < r^2 || (x-s)^2+(y-s)^2 < r^2
        # dielectric inside the circle
        return mu1
    else
        # air outside the circle
        return 1
    end
end

a1 = [1, 0]  # first lattice vector
a2 = [0, 1]  # second lattice vector
#a2 = 0.5*[sqrt(3), 1]

#d1 = 0.005/4  # resolution along first lattice vector
#d2 = 0.005/4  # resolution along second lattice vector

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]

geometry_FDFD = Geometry(epf, muf, a1, a2, d1, d2, TM)
#solver_FDFD = Peacock.FDFD.Solver(geometry_FDFD, [2*d1, 2*d2])
#figure(size(3,4)); plot_band_diagram(solver_FDFD, ks, TM, bands=1:6, dk=0.4, frequency_scale=1/2pi)

Chern_number(epf, muf, a1, a2, d1, d2, TM, Nkx=Nkx, Nky=Nky, bands=bands)
end
#=
geometry_PWE = Geometry(epf, muf, a1, a2, d1, d2)
solver_PWE = Solver(geometry_PWE, 7)
figure(); plot_band_diagram(solver_PWE, ks, TM, bands=1:4, dk=0.4, frequency_scale=1/2pi)
=#