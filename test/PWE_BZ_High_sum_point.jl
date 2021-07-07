using Peacock
using PyPlot, SpecialFunctions, LinearAlgebra
r = 0.2
epf(x, y) = x^2 + y^2 < r^2 ? 9 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.001; d2 = 0.001
geometry = Geometry(epf, muf, a1, a2, d1, d2)
solver = Solver(geometry, 7, 7)
m = 3

kx = [i for i in -m:m]
kx = repeat(kx,1, 2*m+1)
ky = copy(kx)
kx = kx'
K = hcat(kx[:],ky[:])#为何没有*2pi
Kp = K

Kxp = [i for i in Kp[:,1]]
Kxp = repeat(Kxp,1,length(Kp[:,1]))
Kx = copy(Kxp)
Kxp = Kxp'

Kyp = [i for i in Kp[:,2]]
Kyp = repeat(Kyp,1,length(Kp[:,1]))
Ky = copy(Kyp)
Kyp = Kyp'

function geo(kx,ky)
    k=sqrt.(kx.^2 .+ ky.^2)
    eb = 1
    eab=8
    r=0.2
    ep=eab*r*besselj.(1,2*pi*k*r)./k.*cos.(pi*kx)
    ep[k.==0] .= (eb+eab*pi*r^2)
    return ep
end

ep1 = solver.epc
#ep1 = geo(Kx-Kxp, Ky-Kyp)

omegaTM = []
omegaTE = []
for l=0:0.05:1.5;
    #0:0.25:1.5;
    if l<=0.5
        k=[1 0]*l;
    elseif l<=1;
        k=[0.5 l-0.5];
    else
        k=[0.5 0.5]-[1 1]*(l-1);
    end
    
    #find \omega^2
    kk=ones(length(K[:,1]))*k
    #和lecture是等价的，相当于：kx*kx+ky*ky
    omegaTM2, HkTM=eigen(4*pi^2*(kk+K)*((kk+Kp)').*inv(ep1),sortby=real)
    temp = (kk+K)*(kk+Kp)'
    temp = temp[diagind(temp)]
    temp = diagm(temp) 
    omegaTE2, HkTE=eigen(temp,ep1,sortby=real)
    
#    omegaTM = [omegaTM, sqrt.(omegaTM2)]
    push!(omegaTM,sqrt.(Complex.(omegaTM2)))
    push!(omegaTE,sqrt.(Complex.(omegaTE2)))
end
figure()
plot(real.(omegaTM)/2/pi) #这里也没有÷2pi
ylim([0,0.8])