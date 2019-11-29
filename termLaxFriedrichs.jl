include("./genericHam.jl")
include("./SpatialDerivative.jl")
include("./artificialDissipationGLF.jl")

#using .SpatialDerivative

function termLaxFriedrichs(t,y, schemeData)
    data = y
    # Array of type any ~~ similar to cell in MATLAB
    derivL = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivR = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivC = Array{Any, 1}(undef, schemeData.MyGrid.dim)


    # Calculate SpatialDerivative across all the dimensions
    for i = 1:g.dim
        #@time begin
        derivL[i], derivR[i] = upwindFirstFirst(schemeData.MyGrid,data, i)
    #end
        #if i == 3
        #    println(derivL[3])
        #end
        derivC[i] = (derivL[i] + derivR[i])/2
    end

    ham = genericHam(t, data, derivC, schemeData)
    diss, stepBound = artificialDissipationGFL(t, data, derivL, derivR, schemeData)

    if schemeData.uMode == "min"
        ham = -ham
    end

    # (unstable) analytic Ham - disspative stabilization
    delta = ham - diss

    ydot = -delta

    return ydot, stepBound

    #return Ham
end

function odeCFL1(t_interval, y0, schemeData)
    # These parameters can be easily changed
    factorCFL = 0.8
    singleStep = "on"

    y = y0
    t = t_interval[1]
    # Does not matter the time here?
    ydot, stepBound = termLaxFriedrichs(t_interval[1], y0, schemeData)
    deltaT = min(factorCFL* stepBound, t_interval[2] -  t_interval[1])
    t = t + deltaT
    y = y + deltaT.*ydot
    #println(deltaT)
    return t,y
end
