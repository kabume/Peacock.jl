using SparseArrays
"""
    as_to_bs(a1, a2)

Calculate reciprocal lattice vectors, `b1` and `b2`, from the real space lattice
vectors, `a1` and `a2`.
"""
function as_to_bs(a1, a2)
    b1 = 2pi * [+a2[2],-a2[1]] / (+a1[1]*a2[2]-a1[2]*a2[1])
    b2 = 2pi * [-a1[2],+a1[1]] / (-a2[1]*a1[2]+a2[2]*a1[1])
    return b1, b2
end


"""
    bs_to_as(b1, b2)

Calculate real space lattice vectors, `a1` and `a2`, from reciprocal lattice
vectors, `b1` and `b2`.

This is actually the same as `as_to_bs`, but I think having both functions makes
the intention of the code more obvious.
"""
function bs_to_as(b1, b2)
    return as_to_bs(b1, b2)
end

"""
    DEX, DEY, DHX, DHY = diff_yee2(NGRID,RES,BC,kinc=nothing)
=================
- `NGRID`: [Nx Ny] grid size.
- `RES`: [dx dy] grid resolution of the 1X grid.
- `BC`: [xbc ybc] boundary conditions, `-2`: periodic (requires kinc), `0`: Dirichlet.
- `kinc`: [kx ky] incident wave vector. This argument is only needed for periodic boundaries.
Position on the grid: m = (ny - 1)*Nx + nx
"""
function diff_yee2(NGRID::Matrix{Float64}, RES::Matrix{Float64}, BC::Matrix{Int64}, kinc=nothing)
    if kinc != nothing
        kinc0 = kinc
    end
    p(nx,ny,NX)::Int64 = nx+NX*(ny-1)
    dx = RES[1]; dy = RES[2]
    NX = NGRID[1]; NY = NGRID[2]; N = NX*NY;
    Lamx = NX * dx
    Lamy = NY * dy
    N = convert(Int64,N);NX = convert(Int64,NX);NY = convert(Int64,NY)
    if NX == 1
        if BC[1] == 0
            DEX = spzeros(N, N)
        else
            DEX = spdiagm(0 => im*kinc0[2]*ones(N))
        end
    else
        DEX = spdiagm(0 => -1 ./ dx .* ones(ComplexF64, N), 1 => 1 ./ dx .* ones(ComplexF64, N - 1))
        for ny=1:(NY-1)
            DEX[p(NX,ny,NX),p(NX+1,ny,NX)] = 0  # Dirichlet BC
        end
        if BC[1]==-2
            for ny=1:NY
                DEX[p(NX,ny,NX),p(1,ny,NX)]=exp(im*kinc0[2]*Lamx)/dx # periodic BC
            end

        end
    end

    if NY==1
        if BC[2] == 0
            DEY = spzeros(N,N)
        else
            DEY = spdiagm(0 => im*kinc0[1]*ones(N))
        end
    else
        DEY = spdiagm(0 => -1 ./dy.*ones(N), NX => 1 ./dy.*ones(N - NX)) .+ 0*im
        if BC[2]==-2
            for nx=1:NX
                DEY[p(nx,NY,NX),p(nx,1,NX)]=exp(-im*kinc0[1]*Lamy)/dy # periodic BC
            end
        end
    end

    DHX=-DEX'
    DHY=-DEY'
    return DEX,DEY,DHX,DHY
end


# """
#     DiagonalMatrix(diag::AbstractVector{ComplexF64})

# A sparse diagonal matrix that can be used in left division (D \\ X)
# """
# struct DiagonalMatrix <: AbstractMatrix{ComplexF64}
#     diag::AbstractVector{ComplexF64}
# end
# Base.size(A::DiagonalMatrix) = (length(A.diag), length(A.diag))
# Base.getindex(A::DiagonalMatrix, I::Vararg{Int,2}) = I[1]==I[2] ? A.diag[I[1]] : 0
