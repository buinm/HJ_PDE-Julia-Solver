include("./DubinsCar.jl")

using .DubinsCar
function genericPartial(t, data, derivMin, derivMax, schemeData, dim)
    g = schemeData.MyGrid
    dynSys = schemeData.obj

    # Calculate upper and lower bound control
    uU = DubinsCarOptCtrl(schemeData.obj, derivMax, schemeData.uMode)
    uL = DubinsCarOptCtrl(schemeData.obj, derivMin, schemeData.uMode)

    # Calculate disturbance
    dU = DubinsCarOptDstb(schemeData.obj, derivMax, schemeData.dMode)
    dL = DubinsCarOptDstb(schemeData.obj, derivMin, schemeData.dMode)

    # Compute alpha
    dxUU = DunbinsCarDynamics(schemeData.obj, schemeData.MyGrid, uU, dU)
    dxLL = DunbinsCarDynamics(schemeData.obj, schemeData.MyGrid, uL, dL)
    dxLU = DunbinsCarDynamics(schemeData.obj, schemeData.MyGrid, uL, dU)
    dxUL = DunbinsCarDynamics(schemeData.obj, schemeData.MyGrid, uU, dL)


    alpha = max.(abs.(dxUU[dim]), abs.(dxUL[dim]), abs.(dxLL[dim]), abs.(dxLU[dim]))
    return alpha
end
