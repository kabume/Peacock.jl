using Peacock, PyPlot

th1 = 0;
th2 = 60;

xc = 1/4; yc = 1/4/sqrt(3)
#xc = -1/4; yc = 1/4/sqrt(3)
mu1=[0.8736 -0.6671im 0;0.6671im 0.8736 0;0 0 0.8736]
r1 = 0.2
r2 = 0.2
epf(x,y) = (x-xc)^2 + (y-yc)^2 <r1^2 || (x+xc)^2 + (y+yc)^2 < r2^2 ? 15 : 1
muf(x,y) = (x-xc)^2 + (y-yc)^2 <r1^2 || (x+xc)^2 + (y+yc)^2 < r2^2 ? mu1 : 1
d1, d2 = 0.005, 0.005

geometry_FDFD = Geometry(epf, muf, th1, th2, d1, d2, TM)
plot(geometry_FDFD)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
M = BrillouinZoneCoordinate(  0, 1/2, "M")
K = BrillouinZoneCoordinate(-1/3, 1/3, "K")
ks = [G, M, K, G]


solver_FDFD = Peacock.FDFD.Solver(geometry_FDFD, [2*d1, 2*d2])
figure()
@time plot_band_diagram(solver_FDFD, ks, TM, color="blue", bands=1:3, dk=0.2, frequency_scale=1/2pi)