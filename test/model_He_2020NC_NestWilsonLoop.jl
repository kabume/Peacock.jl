using Peacock, PyPlot, LinearAlgebra
orthonormalise = Peacock.orthonormalise
get_k = Peacock.get_k
sample_path = Peacock.sample_path

a = 1; a1 = [1, 0]; a2 = [0, 1]
w = 0.17 # 柱子边长
d = w/2 # w/2 中间， a/2 - w/2 四周 a/4 中心

eps0 = 1; mu0 = 1
mu1 = [14 12.4im 0; -12.4im 14 0; 0 0 1]
eps1 = 15

function epf(x, y)
    if abs(x - d) < w/2 && abs(y - d) < w/2 || abs(x + d) < w/2 && abs(y - d) < w/2 || abs(x - d) < w/2 && abs(y + d) < w/2 || abs(x + d) < w/2 && abs(y + d) < w/2
        return eps1
    else
        return eps0
    end
end

function muf(x, y)
    if abs(x - d) < w/2 && abs(y - d) < w/2 || abs(x + d) < w/2 && abs(y - d) < w/2 || abs(x - d) < w/2 && abs(y + d) < w/2 || abs(x + d) < w/2 && abs(y + d) < w/2
        return mu1
    else
        return mu0
    end
end

#epf(x, y) = x^2 + y^2 < 0.2^2 ? 15 : 1
#muf(x, y) = x^2 + y^2 < 0.2^2 ? 1 : 1

d1 = d2 = 0.005
geometry = Geometry(epf, muf, a1, a2, d1, d2, TM)
solver = Peacock.FDFD.Solver(geometry, 2*[d1, d2])

kx = -pi:0.1:pi
ky = -pi:0.1:pi
modes = Peacock.FDFD.solve(solver,[0,0],TM,bands=1:1)
weighting = modes[1].weighting
data = zeros(ComplexF64, length(kx),length(ky)+1,2,length(weighting))
i = 1
for kx in kx
    j = 1
    for ky in ky
        for n = 1:2
            modes = Peacock.FDFD.solve(solver,[kx,ky],TM,bands=1:2)
            data0 = get_field_FDFD(modes[n].data, modes[n].basis,k0=modes[n].k0)
            data0 = data0[:] / sqrt(abs(dot(data0[:], weighting.*data0[:])))
            data[i,j,n,:] = data0[:]
        end
        j = j + 1
    end
    data[i,end,1,:] = data[i,1,1,:]
    data[i,end,2,:] = data[i,1,2,:]
    i = i + 1
end

vy = zeros(ComplexF64, length(kx),length(ky),2)
W = zeros(ComplexF64, 2,2)
states = zeros(ComplexF64, length(kx),length(ky),2,2)
for i = 1:length(kx)
    for j1 = 1:length(ky)
        F = I
        for j2 = 0:length(ky)-1
            if j1 + j2 <= length(ky)
                j = j1+j2
            else
                j = j1 + j2 - length(ky)
            end
            for k = 1:2
                for l = 1:2
                    W[k,l] = dot(data[i,j,k,:],weighting.*data[i,j+1,l,:])
                end
            end
            F = F * W
        end
        vals,state = eigen(F,sortby=angle)
        vy[i,j1,:] = angle.(vals)
        states[i,j1,:,:] = state 
    end
end

wx = zeros(ComplexF64, length(kx)+1,length(ky),length(weighting))
for j = 1:length(ky) #ky
    for i = 1:length(kx) #kx
        wx[i,j,:] = states[i,j,1,1] * data[i,j,1,:] + states[i,j,2,1] * data[i,j,2,:]
        wx[i,j,:] = wx[i,j,:]/sqrt(abs(dot(wx[i,j,:],weighting.*wx[i,j,:])))
    end
    wx[end,j,:] = wx[1,j,:]
end

px = zeros(length(ky))
for j = 1:length(ky) #ky
    W = 1
    for i = 1:length(kx) #kx
        W = W*dot(wx[i,j,:],weighting.*wx[i+1,j,:])
    end
#    vals = eigvals(W)
    px[j] = angle.(W)
#    px[j] = imag(log(W/2/pi))
end
plot(px,"*")
ylim(-pi,pi)
#figure()
#plot(vy,"*")