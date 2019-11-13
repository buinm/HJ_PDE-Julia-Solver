include("./DubinsCar.jl")
include("./grid.jl")

using .Grid
using .DubinsCar

# For now, assume dynSys to be Dunbins Car
# And uMode = 'min'
# with no disturbance yet
function genericHam(t, data, deriv, grid, obj) # Note that deriv is a cell of values
    uMode = "min"
    uOpt = DubinsCarOptCtrl(obj,deriv, uMode)
    # Skip below for now as we need to
    #d = DubinsCarOptDstb()
    dOpt = 0
    hamValue = 0

    # Plug optimal control into dynamics to compute Hamiltonian
    dx = DunbinsCarDynamics(obj, grid, uOpt, dOpt)

    # Now compute the Hamiltonian
    for i = 1:obj.nx
        hamValue = hamValue .+ deriv[i].*dx[i]
    end

    return hamValue
end
