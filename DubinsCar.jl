
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
    DubinsCar(x, wmax, speed,dMax,dims) = new(x, wmax, speed, dMax, dims)
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

function makeDubinsCar(DubinCar)

end
a = DubinsCar(5,6,7,8,10)
