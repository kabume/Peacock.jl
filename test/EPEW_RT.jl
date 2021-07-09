using Peacock
using PyPlot 
p = TM
r = 0.2
epf(x, y) = x^2 + y^2 < r^2 ? 9 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.001; d2 = 0.001
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, 15, 15)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
figure(); plot_band_diagram(solver, ks, p, bands=1:7, frequency_scale=1/2pi)
grid()

omega0 = 0.005:0.005:1
ky = 0
ii = 1
RT = zeros(ComplexF64, length(omega0))
TT = zeros(ComplexF64, length(omega0))
for omega in omega0
    kT, dT, VT = solver_EPEW(solver, ky, omega, p)
    RT[ii], TT[ii] = abcd(kT, dT, VT, omega, ky)
    ii = ii + 1
end
figure();plot(omega0,abs.(RT))
#plot(omega0, real((1 .+ RT)./(1 .- RT)))
xlabel("ω(2πc/a)")
ylabel("R")
#ylim([-20,20])
grid()

function print_k(omega0)
    for omega = omega0
    kym = (-1:0.05:1)*0.5
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
    M = size(kTD)[2]
    figure()
    plot(kTD[:],(kym*ones(1,M))[:],"r.")
    plot(imag(kTDi[:]),(kym*ones(1,M))[:],"b.")
    xlabel(L"k_x(2\pi/a)")
    ylabel(L"k_y(2\pi/a)")
    title("ω=$omega(2πc/a)")
    legend([L"Re(k_{x})",L"Im(k_x)"])
    grid()
    xlim([-0.6,0.6])
    ylim([-0.6,0.6])
    axhline(0.5,xmin=0.1/1.2,xmax=1.1/1.2,ls="-",color="black",linewidth=0.5)
    axhline(-0.5,xmin=0.1/1.2,xmax=1.1/1.2,ls="-",color="black",linewidth=0.5)
    axvline(0.5,ymin=0.1/1.2,ymax=1.1/1.2,ls="-",color="black",linewidth=0.5)
    axvline(-0.5,ymin=0.1/1.2,ymax=1.1/1.2,ls="-",color="black",linewidth=0.5)
    savefig("r_0p2_eps_9_omega_$omega.svg")
    close()
    end
end
