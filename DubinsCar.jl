include("./grid.jl")

using .Grid
export DunbinsCarDynamics, DubinsCarOptCtrl, DubinsCarOptDstb
mutable struct properties
    xhist
    yhist
end
mutable struct DubinsCar
    x
    wMax
    speed
    dMax
    dims
    DubinsCar(x, wMax, speed) = new(x, wMax, speed)
    # Below variables are for book keeping and plotting
    nx          # Number of state dimensions
    nu          # Number of control inputs
    nd          # Number of disturbance dimensions

    x_state           # State
    u           # Recent control signal

    xhist       # History of state
    uhist       # History of control

    pdim        # position dimensions
    vdim        # velocity dimensions
    hdim        # heading dimensions

    # Figure handle
    hpxpy           # Position
    hpxpyhist       # Position history
    hvxvy           # Velocity
    hvxvyhist       # Velocity history
    # Position velocity (so far only used in DoubleInt)
    hpv
    hpvhist
    # Data (any data that one may want to store for convenience)
    data

end

function makeDubinsCar(x,wMax, speed, )
    # Some basic vehicle properties

    obj.dims = [1;2;3]
    obj.pdim = [1;2]
    obj.hdim = [3]
    obj.nx = 3
    obj.nu = 1
    obj.nd = 3

    obj.xhist = obj.x

end

function DunbinsCarDynamics(obj, grid, u, d)]
    if d == 0 # not using disturbance for now
        d = [0;0;0]
    end
    data = zeros(tuple(grid.pts_each_dim ...,))

    # dx is a cell here
    dx = Array{Any, 1}(undef, 3)
    # Broadcasting using data dimensions
    dx[1] = data .+ (obj.speed .* cos.(grid.vs[1]) .+ d[1])
    dx[2] = data .+ (obj.speed .* sin.(grid.vs[2]) .+ d[2])
    dx[3] = u .+ d[3]

    return dx
end

function DubinsCarOptCtrl(obj, deriv, uMode)
    if uMode == "max"
        uOpt = (deriv[3] .>= 0).*(obj.wMax) + (deriv[3] .< 0).*(-obj.wMax)
    elseif uMode == "min"
        uOpt = (deriv[3] .>= 0).*(-obj.wMax) + (deriv[3] .< 0).*(obj.wMax)
    end

    return uOpt
end

function DubinsCarOptDstb()

end
