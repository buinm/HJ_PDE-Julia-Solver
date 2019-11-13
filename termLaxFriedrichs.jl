
include("./grid.jl")
include("./BasicShapes.jl")
include("./BoundaryCondition.jl")
include("./DubinsCar.jl")
include("./SpatialDerivative.jl")

using .Grid
using .BasicShapes
using .BoundaryCondition
using .DubinsCar
using .SpatialDerivative
# SchemeData variable is unnecesssary  here as all the functions we are using
# only 1 set of utility functions

function termLaxFriedrichs(t,y, g, obj)
    data = reshape(y, tuple(g.pts_each_dim...,)) # kind of redundant?

    # Array of type any ~~ similar to cell in MATLAB
    derivL = Array{Any, 1}(undef, g.dim)
    derivR = Array{Any, 1}(undef, g.dim)
    derivC = Array{Any, 1}(undef, g.dim)


    # Calculate SpatialDerivative across all the dimensions
    for i = 1:g.dim
        derivL[i], derivR[i] = upwindFirstFirst(g,data, i)
        derivC[i] = (derivL[i] + derivR[i])/2
    end

    # TODO: implement these functions :)
    Ham = genericHam(t, data, derivC, g, obj)
    diss, stepBound = artificialDissipationGFL(t, data, derivL, derivR, g, obj)

    # (unstable) analytic Ham - disspative stabilization
    delta = ham - diss

    ydot = -delta

    return ydot, stepBound


end
