#module DubinsCar

include("./grid.jl")
include("./BasicShapes.jl")
include("./SpatialDerivative.jl")

#using .Grid
#using .BasicShapes
#using .SpatialDerivative

export dubinsCar, makeDubinsCar, DunbinsCarDynamics, DubinsCarOptCtrl, DubinsCarOptDstb

mutable struct dubinsCar
    x
    wMax::Float64
    speed::Float64
    dMax
    dubinsCar(x, wMax, speed, dMax) = new(x, wMax, speed, dMax)
    dims
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

function makeDubinsCar(x,wMax, speed, dMax)
    # Some basic vehicle properties
    obj = dubinsCar(x,wMax, speed, dMax)
    obj.dims = [1;2;3]
    obj.pdim = [1;2]
    obj.hdim = [3]
    obj.nx = 3
    obj.nu = 1
    obj.nd = 3

    obj.xhist = obj.x

    return obj
end

function DunbinsCarDynamics(obj, grid, u, d)
    #data = zeros(tuple(grid.pts_each_dim ...,))

    # dx is a cell here
    dx = Array{Any, 1}(undef, 3)
    # Broadcasting using data dimensions
    dx[1] = (obj.speed .* cos.(grid.vs[3]) .+ d[1])
    dx[2] = (obj.speed .* sin.(grid.vs[3]) .+ d[2])
    #println(size(dx[2]))
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

function DubinsCarOptDstb(obj, deriv, dMode)
    dOpt = Array{Any, 1}(undef, obj.nd)
    if dMode == "max"
        for i=1:3
            dOpt[i] = (deriv[i] .>= 0).*(obj.dMax[i]) + (deriv[i] .< 0).*(-obj.dMax[i])
        end
    elseif dMode == "min"
        for i = 1:3
            dOpt[i] = (deriv[i] .>= 0).*(-obj.dMax[i]) + (deriv[i] .< 0).*(obj.dMax[i])
        end
    end
    return dOpt
end

#end
