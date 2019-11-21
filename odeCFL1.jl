include("./termLaxFriedrichs.jl")

function odeCFL1(t_interval, y0, schemeData)
    # These parameters can be easily changed
    factorCFL = 0.8
    singleStep = "on"

    y = y0

    # Does not matter the time here?
    ydot, stepBound = termLaxFriedrichs(t_interval[1], y0, schemeData)
    deltaT = min(factorCFL* stepBound, t_interval[2] -  t_interval[1])

    t = t + deltaT
    y = y + deltaT.*ydot
    return t,y
end
