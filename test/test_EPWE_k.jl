using LinearAlgebra: real
using PyPlot, SpecialFunctions, LinearAlgebra
using Peacock
polarisation = TM
eye(M) = Matrix{ComplexF64}(I, M, M)
diag_R(mat) = mat[diagind(mat)]

function mask(M)
    mat = zeros(M,M)
    for i = -(M-1):M-1
        mat = mat + diagm(i=>(-1)^(i)*ones(M-abs(i)))
    end
    return mat
end

r = 0.3
epf(x, y) = x^2 + y^2 < r^2 ? 9 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.001; d2 = 0.001
#geometry = Geometry(epf,muf,a1,a2,d1,d2)
#solver = Solver(geometry,11,11) #Solver里的Kx，Ky和下面的反过来了，相同的时候还是有问题，应该是a1，a2把方向给固定了

m = 5
M = (2*m+1)^2

kx = [i for i in -m:m]
kx = repeat(kx,1, 2*m+1)
ky = kx'
K = hcat(kx[:],ky[:])#为何没有*2pi
Kp = K
#K = hcat(diag_R(solver.Ky),diag_R(solver.Kx))/2/pi
#Kp = K


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

#ep0 = solver.epc #不一样的原因是Peacock里Kx-Kxp对应的是qs，Ky-Kyp对应的是ps
#ep2 = geo(Ky-Kyp,Kx-Kxp) 
epc = geo(Kx-Kxp,Ky-Kyp) 
#epc = ep0.*mask(size(ep0)[1])
muc = 1

ii = 1
kym = (-1:0.01:1)*0.5
for omega = 0.4
kTD = []
kTDi = []
kTA = []
for ky in kym

    LHS11 = zeros(ComplexF64, M, M)
    LHS12 = eye(M)
    Kp = K + ones(M) * [0 ky]

    if polarisation == TE
        LHS21 = epc * (omega^2*muc*eye(M) .- inv(epc).*(Kp*transpose(Kp)))
        LHS22 = -epc*(inv(epc).*(K[:,1]*ones(1,M)+ones(M)*K[:,1]'))
        LHS = [LHS11 LHS12; LHS21 LHS22]
        RHS = epc
        label = L"H_z"

    elseif polarisation == TM
        LHS21 = (omega^2*epc - diagm(K[:,1].^2 + (ky .+ K[:,2]).^2))
        LHS22 = diagm(-2*K[:,1])
        LHS = [LHS11 LHS12; LHS21 LHS22]
        RHS = muc

        label = L"E_z"
    end
    kT, VT = eigen(LHS,sortby=x -> -abs(x)) #TODO注意这里的排序可能不同, 排序对结果不影响
    #kT, VT = eigen(LHS)
    VT = VT[1:M, :]
    #归一化，不用归一化也行
    VT = VT * diagm(vec(1 ./sqrt.(sum(abs.(VT).^2,dims=1))))
    maxVT,idVT=findmax(abs.(VT), dims=1)
    VT=VT*diagm(vec(maxVT./VT[idVT]))

    dT=1/omega*diag_R(VT'*inv(RHS)*(diagm(K[:,1])*VT+VT*diagm(kT)))

    i0= (abs.(real.(dT)) .> abs.(imag.(dT))).*((real.(dT)+imag.(dT)) .< 0) + (abs.(real.(dT)) .< abs.(imag.(dT))).*(imag.(kT).< 0)
    i0 = Bool.(i0)

    i3 = []
    for i4 = 0:1
        i0 = .!i0
        i1= 1:2*M
        i1 = i1[i0]
        i2 = sortperm(abs.(real.(kT[i0])))
        buf = sort(abs.(real.(kT[i0])))
        i1=i1[i2]
        i2 = (abs.([imag.(kT[i1]); Inf] - [Inf; imag.(kT[i1])] ) .> 1e-10) .| (abs.([imag.(kT[i1]); Inf]) .<= 1e-10) .| (abs.(abs.([real.(kT[i1]); Inf]) .- abs.([Inf; real.(kT[i1])])) .> 1e-10) .| (abs.(abs.([real.(kT[i1]); Inf]) .+ abs.([Inf; real.(kT[i1])]) .- 1) .>= 0.1)
        i1=i1[i2[1:(length(i2)-1)]]
        i3=[i3; i1[1:2m+1]]
    end
    kT=kT[i3]
    dT=dT[i3]
    VT=VT[:,i3]

    buffer=copy(kT)
    i1 = abs.(imag.(buffer)) .< 1e-6
    buffer1=real.(buffer)
    buffer1[.!i1] .= NaN #剔除所有虚部绝对值大于1e-5的值
    if length(buffer1) < length(kT)
        buffer1[length(buffer1)+1:lenm] .= NaN
    end
#    push!(kTD, buffer1')
	kTD = [kTD; buffer1[:]] #所有纯实数解
    buffer[i1] .= NaN
    if length(buffer) < length(kT)
        buffer[length(buffer)+1:lenm] .= NaN
    end
#    push!(kTDi, buffer')
	kTDi = [kTDi; buffer[:]] #所有复数解
    kTA = [kTA; kT[:]]
end
kTD = reshape(kTD,Int(length(kTD)/length(kym)),length(kym))
kTDi = reshape(kTDi,Int(length(kTDi)/length(kym)),length(kym))
kTA = reshape(kTA,Int(length(kTA)/length(kym)),length(kym))
kTD = transpose(kTD)
kTDi = transpose(kTDi)
kTA = transpose(kTA)
figure()
plot(real(kTA[:]),(kym*ones(1,2*(2m+1)))[:],".")
plot(imag(kTA[:]),(kym*ones(1,2*(2m+1)))[:],".")
xlabel(L"k_x(2\pi/a)")
ylabel(L"k_y(2\pi/a)")
title("ω=$omega(2πc/a)")
legend([L"Re(k_{x})",L"Im(k_x)"])
grid()
xlim([-1,1])
savefig("r_0p3_eps_9_omega_$omega.svg")
close()
end
#kTA
#plot(kTD,".")