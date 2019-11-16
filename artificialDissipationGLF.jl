include("./genericPartial.jl")

function artificialDissipationGFL(t, data, derivL, derivR, schemeData)

    derivMin = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivMax = Array{Any, 1}(undef, schemeData.MyGrid.dim)
    derivDiff = Array{Any, 1}(undef, schemeData.MyGrid.dim)

    for i = 1:schemeData.MyGrid.dim
        derivMinL   = minimum(derivL[i])
        derivMinR   = minimum(derivR[i])
        derivMin[i] = min(derivMinL, derivMinR)

        derivMaxL   = maximum(derivL[i])
        derivMaxR   = maximum(derivR[i])
        derivMax[i] = max(derivMaxL, derivMaxR)

        derivDiff[i] = derivR[i] - derivL[i]
    end

    diss = 0
    stepBoundInv = 0
    # This is used for Broadcasting grid.vs
    for i = 1:schemeData.MyGrid.dim
        alpha = genericPartial(t, data, derivMin, derivMax, schemeData, i)
        diss = diss .+ (0.5* derivDiff[i].*alpha)
        stepBoundInv = stepBoundInv + maximum(alpha)/schemeData.MyGrid.dx[i]
    end

    stepBound = 1/ stepBoundInv
    return diss, stepBound

end
