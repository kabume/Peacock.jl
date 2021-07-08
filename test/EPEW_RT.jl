using Peacock
using PyPlot 
p = TM
r = 0.15
epf(x, y) = x^2 + y^2 < r^2 ? 11.43 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.01; d2 = 0.01
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, 11, 11)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
figure(); plot_band_diagram(solver, ks, p, bands=1:3, frequency_scale=1/2pi)
grid()

omega0 = 0.1:0.01:0.5
ky = 0
ii = 1
RT = zeros(length(omega0))
TT = zeros(length(omega0))
for omega in omega0
    kT, dT, VT = solver_EPEW(solver, ky, omega, p)
    RT[ii], TT[ii] = abcd(kT, dT, VT, omega, ky)
    ii = ii + 1
end
figure();plot(omega0,RT)
xlabel("ω(2πc/a)")
ylabel("|r|")
grid()

kym = (-1:0.1:1)*0.5
omega = 0.4
kTD = []
kTDi = []
kTA = []
for ky in kym
    kT, dT, VT = solver_EPEW(solver, ky, omega, p)
	buffer=kT
    i1 = abs.(imag.(buffer)) .< 1e-5
    buffer1=real.(buffer)
    buffer1[.!i1] .= NaN
    if length(buffer1) < length(kT)
        buffer1[length(buffer1)+1:lenm] .= NaN
    end
	kTD = [kTD; buffer1[:]]
    buffer[i1] .= NaN
    if length(buffer) < length(kT)
        buffer[length(buffer)+1:lenm] .= NaN
    end
	kTDi = [kTDi; buffer[:]]
	kTA = [kTA; kT[:]]
end
kTD = reshape(kTD,Int(length(kTD)/length(kym)),length(kym))
kTDi = reshape(kTDi,Int(length(kTDi)/length(kym)),length(kym))
kTA = reshape(kTA,Int(length(kTA)/length(kym)),length(kym))
kTD = transpose(kTD)
kTDi = transpose(kTDi)
kTA = transpose(kTA)
figure()
plot(kTD[:],(kym*ones(1,2*(2m+1)))[:],"r.")
plot(imag(kTDi[:]),(kym*ones(1,2*(2m+1)))[:],"b.")
xlabel(L"k_x(2\pi/a)")
ylabel(L"k_y(2\pi/a)")
title("ω=$omega(2πc/a)")
legend([L"Re(k_{x})",L"Im(k_x)"])
grid()
xlim([-1,1])
#savefig("r_0p3_eps_9_omega_$omega.svg")
#close()