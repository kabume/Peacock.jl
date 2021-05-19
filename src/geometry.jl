@enum Polarisation TE TM

struct Mu
    muxx::Matrix{ComplexF64}
    muyy::Matrix{ComplexF64}
    muzz::Matrix{ComplexF64}
    muxy::Matrix{ComplexF64}
    muyx::Matrix{ComplexF64}
end
struct Eps
    epsxx::Matrix{ComplexF64}
    epsyy::Matrix{ComplexF64}
    epszz::Matrix{ComplexF64}
    epsxy::Matrix{ComplexF64}
    epsyx::Matrix{ComplexF64}
end

struct Geometry
    a1::Array{<:Real,1}
    a2::Array{<:Real,1}
    ep::Eps
    mu::Mu
end

function MaterialTensor(mat::Union{AbstractVecOrMat, AbstractFloat, Int})
    if typeof(mat) <: Int || typeof(mat) <: AbstractFloat
        return Matrix(mat*I, 3, 3)
    elseif typeof(mat) <: AbstractVector && size(mat) == (3, )
        return diagm(mat)
    elseif typeof(mat) <: AbstractMatrix && size(mat) == (3, 3)
        return mat
    else
        throw("only support 3X3 Matrix, 3X1 Vector or a number")
    end
end

function Geometry(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real)
    P = ceil(Int, norm(a1)/(d1 - 1e-6))
    Q = ceil(Int, norm(a2)/(d2 - 1e-6))
    ps = range(-0.5, stop=0.5, length=P)[1:end-1] #why end-1?
    qs = range(-0.5, stop=0.5, length=Q)[1:end-1]
    ps = ps .+ step(ps)/2 # ensure that p and q are centered on zero
    qs = qs .+ step(qs)/2
    xs = [p*a1[1]+q*a2[1] for p in ps, q in qs]
    ys = [p*a1[2]+q*a2[2] for p in ps, q in qs]

    ep = epf.(xs, ys)
    mu = muf.(xs, ys)
    epsij(i, j) = [ep[x, y][i, j] for (x, y) in Iterators.product(1:size(ep)[1], 1:size(ep)[2])]
    muij(i, j) = [mu[x, y][i, j] for (x, y) in Iterators.product(1:size(mu)[1], 1:size(mu)[2])]
    epsxx, epsyy, epszz, epsxy, epsyx = epsij(1, 1), epsij(2, 2), epsij(3, 3), epsij(1, 2), epsij(2, 1)
    muxx, muyy, muzz, muxy, muyx = muij(1, 1), muij(2, 2), muij(3, 3), muij(1, 2), muij(2, 1)

    ep = Eps(epsxx, epsyy, epszz, epsxy, epsyx)
    mu = Mu(muxx, muyy, muzz, muxy, muyx)
    return Geometry(a1, a2, ep, mu)
end

function Geometry(epf::Function, muf::Function, a1_deg::Real, a2_deg::Real, d1::Real, d2::Real)
    a1 = [cosd(a1_deg), sind(a1_deg)]
    a2 = [cosd(a2_deg), sind(a2_deg)]
    return Geometry(epf, muf, a1, a2, d1, d2)
end

function Geometry(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real, polarisation::Polarisation)
    if polarisation == TE
        epfTE(x, y) = inv(MaterialTensor(epf(x, y)))
        mufTE(x, y) = MaterialTensor(muf(x, y))
        return Geometry(epfTE, mufTE, a1, a2, d1, d2)
    else
        epfTM(x, y) = MaterialTensor(epf(x, y))
        mufTM(x, y) = inv(MaterialTensor(muf(x, y)))
        return Geometry(epfTM, mufTM, a1, a2, d1, d2)
    end
end