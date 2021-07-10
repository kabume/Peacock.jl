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
    RT[ii], TT[ii], _, _ = abcd(kT, dT, VT, omega, ky)
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


function psi0(x, y, kin, kr, kt, V, r, t)
    psir = []
    for x0 in x
        buf = []
        if x0 <= 0
            for y0 in y
                buf = [buf; exp(im*(kin[1]*x0+kin[2]*y0))+exp(im*(kr[:,1]*x0+kr[:,2]*y0)).'*r] 
            end
        else
            for y0 = y
                buf = [buf; exp(im*(K[:,1]*x0+K[:,2]*y0)).'*V*diagm(exp(im*(kt*x0+kin[2]*y0)))*t]
            end
        end
        psir=[psir buf]
    end
end

function incident(omega, ky0, dky, N, x, y, p)
theta0=asin(ky0/omega)
if abs(ky0+dky) >= omega
    dky = sign(ky0+dky)*omega
else
    dky = ky0+dky
end
theta1 = asin(dky/omega)
theta = theta1 .+ (0:2*N) ./ (N+eps()) .* (theta0-theta1)
k = omega.*[cos.(theta) sin.(theta)]            #The incident wave's k
weight = abs.(abs.(-N:N) .- N .- 1 .- eps()) ./ (N+1+eps())    #The weight of per incident wave

psir=zeros(ComplexF64, length(y),length(x))
m=7
ii = 1
for ky = k[:,2]
     #Solve the TE and TM Mode
    kT, dT, VT = solver_EPEW(solver, ky, omega, p)
    RT, TT, rT, tT = abcd(kT, dT, VT, omega, ky)
    kr=[-sqrt.(Complex.((omega^2 .- ((-7:7) .+ ky).^2))) ((-7:7) .+ ky)]
    psir = psir+psi0(x, y, k[ii,:], kr, kT[1:(2*m+1)], VT[:,1:(2*m+1)], rT, tT)*weight[ii]

    ii = ii + 1
end

x0 = repeat(x, 1, length(y)); x0 = transpose(x0)
y0 = repeat(y, 1, length(x))
figure()
imshow(x0, y0, abs.(psir).^2)
#pcolor(x0,y0,abs(psir).^2);line([0 0],[y(1) y(length(y))]);
end