include("../grid.jl")
include("../DubinsCar.jl")
include("../BasicShapes.jl")
include("./odeCFL1.jl")
#include("../termLaxFriedrichs.jl")

using Parameters
using GR


@with_kw mutable struct SchemeData
    MyGrid::grid
    obj::dubinsCar
    accuracy = "low"
    uMode = "min"
    dMode = "max"

    end
function odeCFL1()
    # Specify time vector
    t0 = 0
    tMax = 1.0
    dt = 0.05
    tau = t0:dt:tMax

    # Give me a grid to hold data
    g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[100,100,100],3)
    # Dynamics system ~~~ Dubins Car
    my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0, [0.0,0.0,0.0])
    #   Target set
    my_data = CyclinderShape(g, 3.0,[0;0;0],1.0)
    # Data Scheme
    myScheme = SchemeData(MyGrid = g, obj=my_car)

    # initialize memory here
    result = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    debug_grad  = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    derivDiff1 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    derivDiff2 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    derivDiff3 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])

    # Main loop to compute value function V(x,t)
    start = 2
    #small = 1e-4
    count = 0
    V_function = copy(my_data)
    for i = start:length(tau)
        tNow = tau[i-1]
        while tNow < tau[i] - 1e-4
            @time begin
            new_result = odeCFL1([tNow, tau[i]], V_function, myScheme, result, derivDiff1, derivDiff2, derivDiff3, debug_grad)
            end
            tNow = new_result[1]
            V_function = new_result[2]
            count = count +1
        end
    end
    println(count)
    #return debug_grad
    return V_function
end
V_function = odeCFL1()
isosurface(V_function, isovalue=0, rotation = 180)
