using Peacock
a1 = [1, 0]; a2 = [1/2, sqrt(3)/2]
#x = 1 / sqrt(3)
#a1 = 2*[x*cosd(-30), x*sind(-30)]
#a2 = 2*[x*cosd(+30), x*sind(+30)]
function eppf(x,y)
    a1 = [1, 0]; a2 = [1/2, sqrt(3)/2]
    u1, u2, u3 = vcat(a1, 0)', vcat(a2, 0)', [0, 0, 1]'
    g=[u1*u1' u1*u2' u1*u3';u2*u1' u2*u2' u2*u3';u3*u1' u3*u2' u3*u3']
    r = 0.2
    xc1, yc1 = 1/4, 1/6
    xc2, yc2 = -1/4, -1/6
    c1 = [1/4/sqrt(3) 1/4 0]*g
    c2 = [-1/4/sqrt(3) -1/4 0]*g
    xc1 = c1[1]; yc1 = c1[2]
    xc2 = c2[1]; yc2 = c2[2]
    if (x-xc1)^2+(y-yc1)^2 < r^2 || (x-xc2)^2+(y-yc2)^2<r^2
        return 2
    else
        return 1
    end

end

function epf(x, y)
    angles=[60n for n in 0:5]
    angles = [30 210]
    ep_bg = 1
    ep_cyl = 15
    Rs = 1/3
    d1s = d2s = 2Rs/3
    if length(Rs) == 1
        Rs = fill(Rs[1], length(angles))
    end
    if length(d1s) == 1
        d1s = fill(d1s[1], length(angles))
    end
    if length(d2s) == 1
        d2s = fill(d2s[1], length(angles))
    end
    @assert length(Rs) == length(d1s) == length(d2s) == length(angles)
    # Triangular lattice sites
    lx, ly = 1, 1*sqrt(3)
    sites = [
        (0, 0),
        (-lx, 0),
        (+lx, 0),
        (-lx/2, -ly/2),
        (+lx/2, -ly/2),
        (+lx/2, +ly/2),
        (-lx/2, +ly/2),
    ]
    # Hexagonal ring of cylinders at each site
    for (x0,y0) in sites
        for (R,d1,d2,th) in zip(Rs,d1s,d2s,angles)
            x_, y_ = [cosd(th) sind(th); -sind(th) cosd(th)] * [(x-x0); (y-y0)]
            x_ = (x_-R) / (d1/2)
            y_ = (y_-0) / (d2/2)
            if x_^2 + y_^2 <= 1
                return ep_cyl
            end
        end
    end
    return ep_bg
end

muf(x, y) = 1

G = BrillouinZoneCoordinate(  0,   0, "Γ")
M = BrillouinZoneCoordinate(  0, 1/2, "M")
K = BrillouinZoneCoordinate(1/3, 1/3, "K")
G = [0,0]; M = [0,2*pi/sqrt(3)]; K = 2*pi*[1/3,1/sqrt(3)]
ks = [G, M, K, G]

d1, d2 = 0.01, 0.01

geometry = Geometry(epf, muf, a1, a2, d1, d2)
plot(geometry)

solver_PWE =  Solver(geometry, 7)
figure(figsize=(4,3))
plot_band_diagram(solver_PWE, ks, TM, color="blue"; bands=1:2, dk=0.1, frequency_scale=1/2pi)
#ylim(0,0.8)

G = BrillouinZoneCoordinate(  0,   0, "Γ")
M = BrillouinZoneCoordinate(  0, 1/sqrt(3), "M")
K = BrillouinZoneCoordinate(1/3, sqrt(3)/3, "K")
ks = [G, M, K, G]