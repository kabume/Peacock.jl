using Peacock
using PyPlot, LinearAlgebra
p = TM
r = 0.18
cutoff = 11
epf(x, y) = x^2 + y^2 < r^2 ? 11.43 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.01*100/121; d2 = 0.01*100/121
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, cutoff, cutoff)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
X = BrillouinZoneCoordinate(1/2,   0, "X")
M = BrillouinZoneCoordinate(1/2, 1/2, "M")
ks = [G,X,M,G]
figure(); plot_band_diagram(solver, ks, p, bands=1:3, frequency_scale=1/2pi)
grid()

#compute transmission and reflection
omega0 = 0.005:0.005:0.7
ky = 0
ii = 1
RT = zeros(ComplexF64, length(omega0))
TT = zeros(ComplexF64, length(omega0))
for omega in omega0
    kT, dT, VT = solver_EPEW(solver, ky, omega, p)
    RT[ii], TT[ii], _, _ = abcd(kT, dT, VT, omega, ky)
    ii = ii + 1
end
figure();plot(omega0,(abs.(RT)))
#plot(omega0, real((1 .+ RT)./(1 .- RT)))
xlabel("ω(2πc/a)")
ylabel("R")
#ylim([-20,20])
grid()

# compute complex Bloch wavevector with a fixed frequency
function print_k(omega0)
    for omega = omega0
    kym = (-1:0.05:1)*0.5
    kTD = []
    kTDi = []
    kTA = []
    kTDri = []
    for ky in kym
        kT, dT, VT = solver_EPEW(solver, ky, omega, p)
        buffer=copy(kT)
        i1 = abs.(imag.(buffer)) .< 1e-6 #判断是否为纯实数
        #i2 = Bool.(abs.(real.(buffer)) .< 1e-6) .& Bool.(abs.(real.(buffer) .+ 0.5) .< 1e-6) .& Bool.(abs.(real.(buffer) .- 0.5) .< 1e-6)#判断是否为纯虚数
        i2 = abs.(real.(buffer)) .< 1e-6
        buffer1=real.(buffer)
        buffer1[.!i1] .= NaN #不为纯实数的取值为NaN
        kTD = [kTD; buffer1[:]] #所有实数的kT

        buffer2 = imag.(buffer)
        buffer2[.!i2] .= NaN #不为纯虚数的取值为NaN
        kTDi = [kTDi; buffer2[:]] #所有纯虚数的kT

        i3 = .!i1 .& .!i2
        buffer3 = copy(buffer)
        buffer3[.!i3] .= NaN
        kTDri = [kTDri; buffer3[:]]

        # buffer[i1] .= NaN # 纯实数的取值为NaN
        # kTDi = [kTDi; buffer[:]] #所有虚部不为0的kT，实部不一定为0
        kTA = [kTA; kT[:]]
    end
    kTD = reshape(kTD,Int(length(kTD)/length(kym)),length(kym))
    kTDi = reshape(kTDi,Int(length(kTDi)/length(kym)),length(kym))
    kTDri = reshape(kTDri,Int(length(kTDri)/length(kym)),length(kym))
    kTA = reshape(kTA,Int(length(kTA)/length(kym)),length(kym))
    kTD = transpose(kTD)
    kTDi = transpose(kTDi)
    kTDri = transpose(kTDri)
    kTA = transpose(kTA)
    M = size(kTD)[2]
    figure()
    plot(kTD[:],(kym*ones(1,M))[:],"r.")
    plot(kTDi[:],(kym*ones(1,M))[:],"b.")
#    plot(real(kTDri[:]),(kym*ones(1,M))[:],"g.")
#    plot(imag(kTDri[:]),(kym*ones(1,M))[:],"b.")
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

    figure()
    plot(real(kTDri[:]),(kym*ones(1,M))[:],"r.")
    plot(imag(kTDri[:]),(kym*ones(1,M))[:],"b.")
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
#    savefig("r_0p2_eps_9_omega_$omega.svg")
#    close()
    end
end

# draw the field and try to compute Wilson loop
function psi0(x, y, kin, kr, kt, V, r, t)
    psir = []
    for x0 in x
        buf = []
        if x0 <= 0
            for y0 in y
                buf = [buf; exp.(im*(kin[1]*x0+kin[2]*y0)) .+ transpose(exp.(im*(kr[:,1]*x0+kr[:,2]*y0)))*r] 
            end
        else
            for y0 = y
#                buf = [buf; transpose(exp.(im*(K[:,1]*x0 .+ K[:,2]*y0)))*V*diagm(exp.(im*(kt*x0 .+ kin[2]*y0)))*t]
                buf = [buf; transpose(exp.(im*(K[:,1]*x0 .+ K[:,2]*y0)))*V*t] #u function
            end
        end
        psir=[psir; buf]
    end
    return psir
end
K = hcat(diag_R(solver.Ky),diag_R(solver.Kx))/2/pi

P=zeros(ComplexF64,1465,30)
jj = 1
omega = 0.78214949
ky0 = 0; dky = 0.0; N = 0; 
x = 0:0.05:6; y = 0:0.05:6; 

theta0=asin(Complex(ky0/omega))
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
m=Int((cutoff-1)/2)
ii = 1
for ky = k[:,2]'
    ky = Float64.(ky)
     #Solve the TE and TM Mode
    kT, dT, VT = solver_EPEW(solver, ky, omega, p, Norm=false)
#    RT, TT, rT, tT = abcd(kT, dT, VT, omega, ky, background=solver.epc.*mask(size(solver.epc)[1]))
    RT, TT, rT, tT = abcd(kT, dT, VT, omega, ky)
    kr=[-sqrt.(Complex.(omega^2 .- ((-m:m) .+ ky).^2)) ((-m:m) .+ ky)]
    temp = psi0(x, y, k[ii,:], kr, kT[1:(2*m+1)], VT[:,1:(2*m+1)], rT, tT)
    psir = psir .+ reshape(temp, length(y), length(x)) .* weight[ii]
    ii = ii + 1
end
#P[:,jj] = psir[1:10:end]
#交叠积分结果显示 angle(psir_0p782'*psir_0p783)=2.87 angle(psir_0p781'*psir_0p782)=0.169
#jj = jj + 1
#end 

x0 = repeat(x, 1, length(y)); x0 = transpose(x0)
y0 = repeat(y, 1, length(x))
#pcolor(x0,y0,abs.(psir).^2,cmap="hsv")
#画的是E而不是u
figure()
subplot(1,2,1);pcolor(x0,y0,imag.(psir), cmap="coolwarm"); title("imag,ω=$omega(2πc/a)"); colorbar()
subplot(1,2,2);pcolor(x0,y0,real.(psir), cmap="coolwarm"); title("real,ω=$omega(2πc/a)"); colorbar()
#pcolor(x0,y0,abs(psir).^2);line([0 0],[y(1) y(length(y))]);

W = zeros(ComplexF64,29)
for ii = 1:29
    #W[ii] = (P[:,ii] .*geometry.ep.epsxx[1:10:end])'*P[:,ii+1]
    W[ii] = P[:,ii]'*P[:,ii+1]
end
sum(angle.(W))