using Peacock, PyPlot, LinearAlgebra
a = 1; a1 = [1, 0]; a2 = [0, 1]
w = 0.17 # 柱子边长
d = w/2 # w/2 中间， a/2 - w/2 四周 a/4 中心

eps0 = 1; mu0 = 1
mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 1]
eps1 = 15

function epf(x, y)
    if abs(x - d) < w/2 && abs(y - d) < w/2 || abs(x + d) < w/2 && abs(y - d) < w/2 || abs(x - d) < w/2 && abs(y + d) < w/2 || abs(x + d) < w/2 && abs(y + d) < w/2
        return eps1
    else
        return eps0
    end
end

function muf(x, y)
    if abs(x - d) < w/2 && abs(y - d) < w/2 || abs(x + d) < w/2 && abs(y - d) < w/2 || abs(x - d) < w/2 && abs(y + d) < w/2 || abs(x + d) < w/2 && abs(y + d) < w/2
        return mu1
    else
        return mu0
    end
end

#epf(x, y) = x^2 + y^2 < 0.11^2 ? 15 : 1
#muf(x, y) = x^2 + y^2 < 0.11^2 ? mu1 : 1

d1 = d2 = 0.01
geometry = Geometry(epf, muf, a1, a2, d1, d2, TM)
ks = [
    BrillouinZoneCoordinate(-0.5, -0.5, "-X"),
    BrillouinZoneCoordinate(-0.5, 0.0, "Γ"),
    BrillouinZoneCoordinate(-0.5, 0.5, "X")
]
solver = Peacock.FDFD.Solver(geometry, 2*[d1, d2])

dk_inner = nothing
dk_outer = 0.25
bands = 1:2
figure()
plot_wilson_loop_winding(solver, ks, TM, bands, dk_outer=dk_outer, dk_inner=dk_inner,delta_brillouin_zone=BrillouinZoneCoordinate(1,0))