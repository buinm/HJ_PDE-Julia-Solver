include("./grid.jl")
include("./DubinsCar.jl")

using .Grid
using .DubinsCar
using Parameters

@with_kw mutable struct sf
    MyGrid::grid
    obj::dubinsCar
    accuracy = "low"
    uMode = "min"
    dMode = "min"
end

g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0)
myScheme = sf(MyGrid = g, obj=my_car)
