function solver_EPWE(solver::Solver, ky::Number, omega::Number, polarisation::Polarisation; Norm=false)
    epc = solver.epc.*mask(size(solver.epc)[1]) 
    muc = solver.muc.*mask(size(solver.muc)[1])
    global K = hcat(diag_R(solver.Ky),diag_R(solver.Kx))/2/pi
    M = length(K[:,1])
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
    kT, VT = eigen(LHS,sortby=x -> -abs(x))
    VT = VT[1:M, :]
    # Normalization
    if Norm == true
    VT = VT * diagm(vec(1 ./sqrt.(sum(abs.(VT).^2,dims=1))))
    maxVT,idVT=findmax(abs.(VT), dims=1)
    VT=VT*diagm(vec(maxVT./VT[idVT]))
    end
    #direction
    dT=1/omega*diag_R(VT'*inv(RHS)*(diagm(K[:,1])*VT+VT*diagm(kT)))

    # find the solution in 1st BZ
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
        i3=[i3; i1[1:Int(sqrt(M))]]
    end
    kT=kT[i3]
    dT=dT[i3]
    VT=VT[:,i3]
    return kT, dT, VT
end

function abcd(kT::VecOrMat,dT::VecOrMat,VT::VecOrMat,omega::Number,ky::Number;background=1)
    M = size(VT)[1]
    m = Int((sqrt(M)-1)/2)
    Ky=Int(round(ky))
    ky=ky-Ky
    msum=zeros((2*m+1),(2*m+1)^2)
    for i0=1:(2*m+1)
        msum[i0,((i0-1)*(2*m+1)+1):(i0*(2*m+1))] .= 1
    end

    ep = background
    abcd1=[-eye(2*m+1) msum*VT[:,1:(2*m+1)]; diagm(sqrt.(Complex.(omega^2 .- ((-m:m) .+ky).^2))) msum*(inv(ep)*(diagm(K[:,1])*VT[:,1:(2*m+1)]+VT[:,1:(2*m+1)]*diagm(kT[1:(2*m+1)])))]
    input = zeros(ComplexF64, 4*m+2) 
    input[Ky+m+1] = 1
    input[Ky+3*m+2] = sqrt(Complex(omega^2-(ky+Ky)^2))
    output=abcd1\input
    rT=output[1:(2*m+1)]
    tT=output[(2*m+2):(4*m+2)]

    kxi = sqrt(Complex(omega^2 - ky^2))
    if kxi==0
        kxi=eps()
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

#    RT = abs(RT); TT = abs(TT)
    return RT, TT, rT, tT
end