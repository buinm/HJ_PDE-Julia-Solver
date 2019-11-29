function ComputeOptmCtrl(dV_dtheta, schemeData)
    if dV_dtheta < 0
        if myScheme.uMode == min
            u = schemeData.obj.wMax
        else
            u = -schemeData.obj.wMax
        end
    else
        if schemeData.uMode == min
            u = (-schemeData.obj.wMax)
        else
            u = schemeData.obj.wMax
        end
    end

    return u
end

function DubinsCarDynamics(schemeData, theta, u)
    dx = Array{Float32, 1}(undef, 3)
    dx[1] =  schemeData.obj.speed* cos(theta)
    dx[2] =  schemeData.obj.speed* sin(theta)
    dx[3] = u

    return dx
end


function Hamiltonian(dV_dstate, k,myScheme) # dV_dstate is a vector of gradient size 3
    uOpt = ComputeOptmCtrl(dV_dstate[3], myScheme)
    # dx is a vector
    dstate_dt = DubinsCarDynamics(myScheme,myScheme.MyGrid.vs[3][1,1,k] ,uOpt)

    if myScheme.uMode == "max"
        return dV_dstate[1]* dstate_dt[1] + dV_dstate[2]* dstate_dt[2] + dV_dstate[3]* dstate_dt[3]
    else
         return -(dV_dstate[1]* dstate_dt[1] + dV_dstate[2]* dstate_dt[2] + dV_dstate[3]* dstate_dt[3])
    end
end

function addGhostExtrapolate(left, right, mode)
    if mode == "Bottom"
        slope = abs(right - left) * 1 *sign(left)
        extended_point = left + slope
        return extended_point
    elseif mode == "Top"
        slope = abs(right - left) *1 * sign(right)
        extended_point = right + slope
        return extended_point
    end
end

"""function addGhostPeriodic(left, right, mode)
    if mode == "Bottom"
        extended_point =
    else if mode == "Top"

    end
end"""
function spatial_Derivative(grids, index ,data_array, dim)
    dxInv = 1/grids.dx[dim]

    if index != 1 && index != grids.pts_each_dim[dim]
        dV_d_L = (data_array[index] - data_array[index-1])*dxInv
        dV_d_R = (data_array[index + 1] - data_array[index])*dxInv
        return dV_d_L, dV_d_R
    else
        if index == 1
            if dim != grids.pDims
                left_point = addGhostExtrapolate(data_array[index], data_array[index + 1], "Bottom")
            else
                #left_point = addGhostPeriodic(data_array[index], data_array[index + 1], "Bottom")
                left_point = data_array[grids.pts_each_dim[dim]]
            end
            dV_d_L = (data_array[index] - left_point)*dxInv
            dV_d_R = (data_array[index + 1] - data_array[index])*dxInv
            return dV_d_L, dV_d_R
        else # index == grids.pts_each_dim[dim]
            if dim != grids.pDims
                right_point = addGhostExtrapolate(data_array[index-1], data_array[index], "Top")
            else
                #right_point = addGhostPeriodic(data_array[index-1], data_array[index], "Top")
                right_point = data_array[1]
            end
            dV_d_L = (data_array[index] - data_array[index-1])*dxInv
            dV_d_R = (right_point- data_array[index])*dxInv
            return dV_d_L, dV_d_R
        end
    end

end

function odeCFL1(tspan, V_function, myScheme)

    # Calculate the Hamiltonian terms
    #global V_function
    #print(V_function)

    min_deriv_1 = 1e9
    min_deriv_2 = 1e9
    min_deriv_3 = 1e9

    max_deriv_1 = -1e9
    max_deriv_2 = -1e9
    max_deriv_3 = -1e9

    t = tspan[1]

    #result = zeros(tuple(myScheme.MyGrid.pts_each_dim ...,))
    result = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    #derivDiff = Array{Any, 1}(undef, myScheme.MyGrid.dim)
    derivDiff1 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    derivDiff2 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    derivDiff3 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    for i = 1:myScheme.MyGrid.pts_each_dim[1] #x index
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for k = 1:myScheme.MyGrid.pts_each_dim[3] #theta
                # Left and Right spatial derivative
                dV_dx_L, dV_dx_R = spatial_Derivative(myScheme.MyGrid, i, V_function[:, j,k], 1)
                dV_dy_L, dV_dy_R = spatial_Derivative(myScheme.MyGrid, j, V_function[i, :,k], 2)
                dV_dtheta_L, dV_dtheta_R = spatial_Derivative(myScheme.MyGrid, k, V_function[i, j,:], 3)
                # Average Spatial Derivative
                dV_dx_C = (dV_dx_L + dV_dx_R)/2
                dV_dy_C = (dV_dy_L + dV_dy_R)/2
                dV_dtheta_C = (dV_dtheta_L + dV_dtheta_R)/2

                dV_dState= [dV_dx_C dV_dy_C dV_dtheta_C]
                ham = Hamiltonian(dV_dState, k ,myScheme)
                result[i,j,k] =  ham

                # Keep track of smallest and largest gradient for dissipation calculation
                min_deriv_1 = min(dV_dx_L, dV_dx_R, min_deriv_1)
                min_deriv_2 = min(dV_dy_L, dV_dy_R, min_deriv_2)
                min_deriv_3 = min(dV_dtheta_L, dV_dtheta_R, min_deriv_3)

                max_deriv_1 = max(dV_dx_L, dV_dx_R, max_deriv_1)
                max_deriv_2 = max(dV_dy_L, dV_dy_R, max_deriv_2)
                max_deriv_3 = max(dV_dtheta_L, dV_dtheta_R, max_deriv_3)

                # Keep track of gradient difference for dissipation calculation
                derivDiff1[i,j,k] =  dV_dx_R - dV_dx_L
                derivDiff2[i,j,k] =  dV_dy_R - dV_dy_L
                derivDiff3[i,j,k] =  dV_dtheta_R - dV_dtheta_L

            end
        end
    end

    delta_t = tspan[2] - tspan[1]
    t = t + delta_t


    # Calculate Dv_dX
    #dV_dx = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],1,1)
    dV_dx = Array{Float64, 3}(undef, 1,1,myScheme.MyGrid.pts_each_dim[1])
    dV_dy = Array{Float64, 3}(undef, 1,1, myScheme.MyGrid.pts_each_dim[2])

    for i = 1:myScheme.MyGrid.pts_each_dim[3]
        theta = myScheme.MyGrid.vs[3][1,1,i]
        dV_dx[1,1,i] = abs(myScheme.obj.speed* cos(theta))
        dV_dy[1,1,i] = abs(myScheme.obj.speed* sin(theta))
    end

    # Calculate dissipation,
    diss = 0.5.*(derivDiff1.*dV_dx + derivDiff2.*dV_dy + derivDiff3)

    for i = 1:myScheme.MyGrid.pts_each_dim[1] #x index
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for k = 1:myScheme.MyGrid.pts_each_dim[3]
                result[i,j,k] = diss[i,j,k] -result[i,j,k]
            end
        end
    end

    V_function = V_function + result .* delta_t

    return  t, V_function




end
