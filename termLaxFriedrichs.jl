include("./genericHam.jl")
include("./SpatialDerivative.jl")
include("./artificialDissipationGLF.jl")

using .SpatialDerivative

function termLaxFriedrichs(t,y, schemeData)
    data = y
    # Array of type any ~~ similar to cell in MATLAB
    derivL = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivR = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivC = Array{Any, 1}(undef, schemeData.MyGrid.dim)


    # Calculate SpatialDerivative across all the dimensions
    for i = 1:g.dim
        derivL[i], derivR[i] = upwindFirstFirst(schemeData.MyGrid,data, i)
        derivC[i] = (derivL[i] + derivR[i])/2
    end

    # TODO: implement these functions :)
    ham = genericHam(t, data, derivC, schemeData)
    diss, stepBound = artificialDissipationGFL(t, data, derivL, derivR, schemeData)

    # (unstable) analytic Ham - disspative stabilization
    delta = ham - diss

    ydot = -delta

    return ydot, stepBound

    #return Ham
end
