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

r = 0.2
epf(x, y) = x^2 + y^2 < r^2 ? 9 : 1
muf(x,y) = 1
a1 = [1, 0]; a2 = [0, 1]
d1 = 0.001; d2 = 0.001
geometry = Geometry(epf,muf,a1,a2,d1,d2)
solver = Solver(geometry,11,11) #Solver里的Kx，Ky和下面的反过来了，相同的时候还是有问题，应该是a1，a2把方向给固定了

m = 5
M = (2*m+1)^2

kx = [i for i in -m:m]
kx = repeat(kx,1, 2*m+1)
ky = kx'
K = hcat(kx[:],ky[:])#为何没有*2pi
Kp = K
K = hcat(diag_R(solver.Ky),diag_R(solver.Kx))/2/pi
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

ep0 = solver.epc #不一样的原因是Peacock里Kx-Kxp对应的是qs，Ky-Kyp对应的是ps
#ep2 = geo(Ky-Kyp,Kx-Kxp) 
#epc = geo(Kx-Kxp,Ky-Kyp) 
epc = ep0.*mask(size(ep0)[1])
muc = 1

omega0 = 0.7
RA = zeros(length(omega0))
TA = copy(RA)
ii = 1
for omega in omega0

#omega = 0.2
ky = 0
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
#kT, VT = eigen(LHS,sortby=x -> -abs(x)) #TODO注意这里的排序可能不同, 排序对结果不影响
kT, VT = eigen(LHS)
VT = VT[1:M, :]
# 归一化，不用归一化也行
#VT = VT * diagm(vec(1 ./sqrt.(sum(abs.(VT).^2,dims=1))))
#maxVT,idVT=findmax(abs.(VT), dims=1)
#VT=VT*diagm(vec(maxVT./VT[idVT]))

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

Ky=Int(round(ky))
ky=ky-Ky
msum=zeros((2*m+1),(2*m+1)^2)
for i0=1:(2*m+1)
    msum[i0,((i0-1)*(2*m+1)+1):(i0*(2*m+1))] .= 1
end

ep = 1
abcd1=[-eye(2*m+1) msum*VT[:,1:(2*m+1)]; diagm(sqrt.(Complex.(omega^2 .- ((-m:m) .+ky).^2))) msum*(inv(ep)*(diagm(K[:,1])*VT[:,1:(2*m+1)]+VT[:,1:(2*m+1)]*diagm(kT[1:(2*m+1)])))]
input = zeros(4*m+2) 
input[Ky+m+1] = 1
input[Ky+3*m+2] = sqrt(omega^2-(ky+Ky)^2)
output=abcd1\input
rT=output[1:(2*m+1)]
tT=output[(2*m+2):(4*m+2)]

kxi = sqrt(omega^2 - ky^2)
if kxi==0
    kxi=eps
end
kxr=sqrt.(Complex.(omega^2 .- (ky - round(ky) .+ (-m:m)).^2))

i0=abs.(imag.(kxr[1:(2*m+1)])) .<= 1e-15
#    i0 = Bool.[i0]
if length(rT[i0])==0
    RT = 0
else
    RT = rT[i0]'*diagm(kxr[i0])*rT[i0]/kxi;
end

i0 = abs.(imag.(dT[1:(2*m+1)])) .<= abs.(real.(dT[1:(2*m+1)]))
if length(tT[i0])==0
    TT=0
else
    TT=tT[i0]'*diagm(dT[1:2m+1][i0])*tT[i0]/kxi*omega;
end

RT = abs(RT); TT = abs(TT)

RA[ii] = abs(RT); TA[ii] = abs(TT)
global  ii = ii + 1
end
figure()
plot(omega0,RA)
