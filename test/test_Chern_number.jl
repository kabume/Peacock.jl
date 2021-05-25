using Peacock
r = 0.1
# mu1 = 1
mu1 = [14 12.4im 0;-12.4im 14 0;0 0 1]
epf(x, y) = x^2 + y^2 <= r^2 ? 15 : 1
muf(x, y) = x^2 + y^2 <= r^2 ? mu1 : 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.01; d2 = 0.01
polar = TM
@time Chern_number(epf, muf, a1, a2, d1, d2, polar, bands=1:4)


mu1=[0.8736 -0.6671im 0;0.6671im 0.8736 0;0 0 0.8736]
r1 = 0.2
r2 = 0.2
xc = 1/4; yc = 1/4/sqrt(3)
epf(x,y) = (x-xc)^2 + (y-yc)^2 <r1^2 || (x+xc)^2 + (y+yc)^2 < r2^2 ? 15 : 1
muf(x,y) = (x-xc)^2 + (y-yc)^2 <r1^2 || (x+xc)^2 + (y+yc)^2 < r2^2 ? mu1 : 1
d1, d2 = 0.01, 0.01
a1 = [1, 0]; a2 = [1/2, sqrt(3)/2]
@time Chern_number(epf, muf, a1, a2, d1, d2, polar, Nkx = 4, Nky = 4, bands=1:3)