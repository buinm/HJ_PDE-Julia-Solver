function genericHam(t, data, deriv, schemeData) # Note that deriv is a cell of values
    uOpt = DubinsCarOptCtrl(schemeData.obj, deriv, schemeData.uMode)
    dOpt = DubinsCarOptDstb(schemeData.obj, deriv, schemeData.dMode)

    hamValue = 0

    # Plug optimal control into dynamics to compute Hamiltonian
    dx = DunbinsCarDynamics(schemeData.obj, schemeData.MyGrid, uOpt, dOpt)

    # Now compute the Hamiltonian
    for i = 1:schemeData.obj.nx
        hamValue = hamValue .+ deriv[i].*dx[i]
    end

    return hamValue
end
