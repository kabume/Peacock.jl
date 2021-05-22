using Peacock, PyPlot

th1 = 0;
th2 = 60;

xc = 1/4; yc = 1/4/sqrt(3)
#xc = -1/4; yc = 1/4/sqrt(3)
#mu1=[0.8736 -0.6671im 0;0.6671im 0.8736 0;0 0 0.8736]
mu1 = 1
r = 0.2
epf(x,y) = (x-xc)^2 + (y-yc)^2 <r^2 || (x+xc)^2 + (y+yc)^2 < r^2 ? 15 : 1
muf(x,y) = (x-xc)^2 + (y-yc)^2 <r^2 || (x+xc)^2 + (y+yc)^2 < r^2 ? mu1 : 1
d1, d2 = 0.01, 0.01

geometry = Geometry(epf, muf, th1, th2, d1, d2)
#plot(geometry)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
M = BrillouinZoneCoordinate(  0, 1/2, "M")
K = BrillouinZoneCoordinate(-1/3, 1/3, "K")
ks = [G, M, K, G]

solver_PWE =  Solver(geometry, 7)
figure(figsize=(4,3))
@time plot_band_diagram(solver_PWE, ks, TM, color="blue"; bands=1:2, dk=0.1, frequency_scale=1/2pi)

geometry_FDFD = Geometry(epf, muf, th1, th2, d1, d2, TM)
solver_FDFD = Peacock.FDFD.Solver(geometry_FDFD, [2*d1, 2*d2])
@time plot_band_diagram(solver_FDFD, ks, TM, color="red", bands=1:2, dk=0.1, frequency_scale=1/2pi)