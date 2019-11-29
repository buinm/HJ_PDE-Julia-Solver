include("../grid.jl")
include("../DubinsCar.jl")
include("../BasicShapes.jl")
include("./odeCFL1.jl")
#include("../termLaxFriedrichs.jl")

#using .Grid
#using .DubinsCar
#using .BasicShapes
using Parameters
using GR


@with_kw mutable struct SchemeData
    MyGrid::grid
    obj::dubinsCar
    accuracy = "low"
    uMode = "min"
    dMode = "max"
end

# Specify time vector
t0 = 0
tMax = 1.0
dt = 0.05
tau = t0:dt:tMax

# Give me a grid to hold data
g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[100,100,100],3)
# Dynamics system ~~~ Dubins Car
my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0, [0,0,0])
#   Target set
my_data = CyclinderShape(g, 3.0,[0;0;0],1.0)
# Data Scheme
myScheme = SchemeData(MyGrid = g, obj=my_car)

#t, V = odeCFL1([0.5, 0.75], my_data,myScheme)

# Main loop to compute value function V(x,t)

start = 2
small = 1e-4
global count = 0
global V_function = copy(my_data)
for i = start:length(tau)
    tNow = tau[i-1]
    while tNow < tau[i] - 1e-4
        global V_function
        global count
        @time begin
        result = odeCFL1([tNow, tau[i]], V_function, myScheme)
        end
        tNow = result[1]
        #println(tNow)
        V_function = result[2]
        count = count +1
    end
end
println(count)

isosurface(V_function, isovalue=0)