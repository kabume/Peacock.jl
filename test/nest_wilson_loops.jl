function plot_wilson_loop_winding(solver::Peacock.FDFD.Solver, ks, polarisation, bands::AbstractVector{<:Int};
    dk_outer=nothing, dk_inner=nothing, delta_brillouin_zone=BrillouinZoneCoordinate(0,1),
    labels=[], markersize=nothing)

    if labels == []
        labels = [hasproperty(x,:label) ? x.label : "" for x in ks]
    end
    ks = [typeof(x)==BrillouinZoneCoordinate ? get_k(x,solver.basis) : x for x in ks]
    ks_outer, _ = sample_path(ks, dk=dk_outer)
    function get_data(ks_outer) #not efficient
        spaces = []
        for k in ks_outer
            spaces_inner = HilbertSpace_FDFD[]
            delta_k = get_k(delta_brillouin_zone, solver.basis)
            ks_inner = [k, k+delta_k]
            ks_inner, _ = sample_path(ks_inner, dk=dk_inner)
            for k in ks_inner
                modes = Peacock.FDFD.solve(solver, k, polarisation)
                space = HilbertSpace_FDFD(modes[bands])
                push!(spaces_inner, space)
            end
            spaces_inner = [spaces_inner; spaces_inner[1]]
            push!(spaces, spaces_inner)
        end
        return spaces
    end
    plot_band_diagram(my_solve, ks, dk=dk_outer, labels=labels, show_vlines=false, markersize=markersize)
    ylim(-pi,pi)
    yticks((-1:1)*pi, labels=["-π","0","+π"])
    ylabel("Wilson spectrum")
end

function get_data2(kx, ky, weighting; bands=1:2)
    data = zeros(ComplexF64, length(kx),length(ky)+1, length(bands), length(weighting))
    for i in 1:length(kx), j in 1:length(ky), n in 1:length(bands)   
        modes = Peacock.FDFD.solve(solver,[kx[i],ky[j]],TM,bands=bands)
        data0 = get_field_FDFD(modes[n].data, modes[n].basis,k0=modes[n].k0)
        data0 = data0[:] / sqrt(abs(dot(data0[:], weighting.*data0[:])))
        data[i,j,n,:] = data0[:]
        if j == length(ky)
            data[i,end,n,:] = data[i,1,n,:]
        end
    end
    return data
end

@everywhere using Peacock
function get_data(kx, ky, weighting; bands=1:2)
    data = zeros(ComplexF64, length(kx),length(ky)+1, length(bands), length(weighting))
    @sync begin
        @distributed for i in 1:length(kx)
            for j in 1:length(ky)
                for n in 1:length(bands)   
                    modes = Peacock.FDFD.solve(solver,[kx[i],ky[j]],TM,bands=bands)
                    data0 = get_field_FDFD(modes[n].data, modes[n].basis,k0=modes[n].k0)
                    data0 = data0[:] / sqrt(abs(dot(data0[:], weighting.*data0[:])))
                    data[i,j,n,:] = data0[:]
                    if j == length(ky)
                        data[i,end,n,:] = data[i,1,n,:]
                    end
                end
            end
        end
    end
    return data
end

function get_data(kx, ky, weighting; bands=1:2)
    data = zeros(ComplexF64, length(kx),length(ky)+1, length(bands), length(weighting))
    inds = CartesianIndices((1:length(kx),1:length(ky),1:length(bands)))
    @sync begin
        @distributed for ind in inds
            i, j, n = ind.I             
            modes = Peacock.FDFD.solve(solver,[kx[i],ky[j]],TM,bands=bands)
            data0 = get_field_FDFD(modes[n].data, modes[n].basis,k0=modes[n].k0)
            data0 = data0[:] / sqrt(abs(dot(data0[:], weighting.*data0[:])))
            data[i,j,n,:] = data0[:]
            if j == length(ky)
                data[i,end,n,:] = data[i,1,n,:]
       
            end
        end
    end
    return data
end