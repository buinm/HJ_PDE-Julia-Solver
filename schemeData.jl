include("./grid.jl")
include("./DubinsCar.jl")
include("./BasicShapes.jl")
include("./termLaxFriedrichs.jl")

using .Grid
using .DubinsCar
using .BasicShapes
using Parameters

@with_kw mutable struct SchemeData
    MyGrid::grid
    obj::dubinsCar
    accuracy = "low"
    uMode = "min"
    dMode = "min"
end

g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0, [0,0,0])
my_data = CyclinderShape(g, 3.0,[0;0;0],1.0)
myScheme = SchemeData(MyGrid = g, obj=my_car)

ham, burger = termLaxFriedrichs(2, my_data,myScheme)
