using Peacock
function Chern_number(solver::Peacock.FDFD.Solver, polar::Polarisation, Nkx = 4, Nky = 4, bands = 1:4)
    a1 = solver.basis.a1
    a2 = solver.basis.a2
    Px = 1     # lattice constant along x
    Py = 1     # lattice constant along y
    N_max = bands[end]  # maximum band number ***

    Pkx = 2pi/Px/a2[2]
    Pky = 2pi/Py/a2[2]
    deltakx = Pkx/Nkx
    deltaky = Pky/Nky
    NK = Nkx*Nky
    kcx = zeros(NK, 4)
    kcy = zeros(NK, 4)
    fieldTot = zeros(N_max, 1)
    fieldtemp = zeros(ComplexF64, N_max, 1)
    
    for m=1:Nkx
        for n=1:Nky
            
            index = (n-1)*Nkx+m
            kcxg = [(m-1)*deltakx,(m)*deltakx,(m)*deltakx,(m-1)*deltakx]
            kcyg = [(n-1)*deltaky,(n-1)*deltaky,(n)*deltaky,(n)*deltaky]
            kcx[index, :] = kcxg * a2[2] + kcyg * a1[2]
            kcy[index, :] = kcyg * a1[1] - a2[1]*kcxg  # Rotation coordinates

            modes = Peacock.FDFD.solve(solver, [kcx[index, 1], kcy[index, 1]], polar, bands = 1:N_max)
            field1 = [get_field_FDFD(modes[i].data, modes[i].basis, k0=modes[i].k0) for i in 1:N_max]

            modes = Peacock.FDFD.solve(solver, [kcx[index, 2], kcy[index, 2]], polar, bands = 1:N_max)
            field2 = [get_field_FDFD(modes[i].data, modes[i].basis, k0=modes[i].k0) for i in 1:N_max]

            modes = Peacock.FDFD.solve(solver, [kcx[index, 3], kcy[index, 3]], polar, bands = 1:N_max)
            field3 = [get_field_FDFD(modes[i].data, modes[i].basis, k0=modes[i].k0) for i in 1:N_max]

            modes = Peacock.FDFD.solve(solver, [kcx[index, 4], kcy[index, 4]], polar, bands = 1:N_max)
            field4 = [get_field_FDFD(modes[i].data, modes[i].basis, k0=modes[i].k0) for i in 1:N_max]
            
            weighting = transpose(modes[1].weighting)

            for mm = 1:N_max
                temp = field1[mm][:]' .* weighting * field2[mm][:] / abs.(field1[mm][:]' .* weighting * field2[mm][:]) * 
                    field2[mm][:]' .* weighting * field3[mm][:] / abs.(field2[mm][:]' .* weighting * field3[mm][:]) *
                    field3[mm][:]' .* weighting * field4[mm][:] / abs.(field3[mm][:]' .* weighting * field4[mm][:]) *
                    field4[mm][:]' .* weighting * field1[mm][:] / abs.(field4[mm][:]' .* weighting * field1[mm][:])

                temp = field1[mm][:]' .* weighting * field2[mm][:] * field2[mm][:]' .* weighting * field3[mm][:] * field3[mm][:]' .* weighting * field4[mm][:] * field4[mm][:]' .* weighting * field1[mm][:] #不归一化对结果没影响
                fieldtemp[mm] = temp[1]
                fieldTot[mm] = fieldTot[mm] + imag(log(fieldtemp[mm]))
            end
        end
    end
    return round.(fieldTot[:]/2/pi)    
end


function Chern_number(epf::Function, muf::Function, a1::Array{<:Real,1}, a2::Array{<:Real,1}, d1::Real, d2::Real, polar::Polarisation; Nkx = 4, Nky = 4, bands=1:4)
    geometry = Geometry(epf, muf, a1, a2, d1, d2, polar)
    solver = Peacock.FDFD.Solver(geometry, 2*[d1, d2])
    return Chern_number(solver, polar, Nkx=Nkx, Nky=Nky, bands = bands)
end
