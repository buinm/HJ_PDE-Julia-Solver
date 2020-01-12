
# Add point at boundary utility function
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

# Function to calculate spatial derivative
function spatial_Derivative(grids, index ,data_array, dim)
    dxInv = 1/grids.dx[dim]

    if index != 1 && index != grids.pts_each_dim[dim]
        dV_d_L = (data_array[index] - data_array[index-1])*dxInv
        dV_d_R = (data_array[index + 1] - data_array[index])*dxInv
        return dV_d_L, dV_d_R
    end
    if index == 1
        if dim != grids.pDims
            left_point = addGhostExtrapolate(data_array[index], data_array[index + 1], "Bottom")
        else
            left_point = data_array[grids.pts_each_dim[dim]]
        end
            dV_d_L = (data_array[index] - left_point)*dxInv
            dV_d_R = (data_array[index + 1] - data_array[index])*dxInv
            return dV_d_L, dV_d_R
    else # index == grids.pts_each_dim[dim]
        if dim != grids.pDims
            right_point = addGhostExtrapolate(data_array[index-1], data_array[index], "Top")
        else
            right_point = data_array[1]
        end
        dV_d_L = (data_array[index] - data_array[index-1])*dxInv
        dV_d_R = (right_point- data_array[index])*dxInv
        return dV_d_L, dV_d_R
    end
end

function odeCFL1(tspan, V_function, myScheme, result, derivDiff1, derivDiff2, derivDiff3)
    # Variables to keep track of the max and min derivative
    min_deriv_1 = 1e9
    min_deriv_2 = 1e9
    min_deriv_3 = 1e9
    max_deriv_1 = -1e9
    max_deriv_2 = -1e9
    max_deriv_3 = -1e9
    thetas = myScheme.MyGrid.min[3] : myScheme.MyGrid.dx[3] : myScheme.MyGrid.max[3]

    # Calculate the Hamiltonian terms
    for k = 1:myScheme.MyGrid.pts_each_dim[3] #theta index
            theta = thetas[k]
            cos_theta = cos(theta)
            sin_theta = sin(theta)
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for i = 1:myScheme.MyGrid.pts_each_dim[1] # x
                # Left and Right spatial derivative
                @views dV_dx_L, dV_dx_R = spatial_Derivative(myScheme.MyGrid, i, V_function[:, j,k], 1)
                @views dV_dy_L, dV_dy_R = spatial_Derivative(myScheme.MyGrid, j, V_function[i, :,k], 2)
                @views dV_dtheta_L, dV_dtheta_R = spatial_Derivative(myScheme.MyGrid, k, V_function[i, j,:], 3)

                # Average Spatial Derivative
                dV_dx_C = (dV_dx_L + dV_dx_R)/2
                dV_dy_C = (dV_dy_L + dV_dy_R)/2
                dV_dtheta_C = (dV_dtheta_L + dV_dtheta_R)/2

                # Calculate optimal control input
                uOpt = myScheme.obj.wMax # Default value
                if dV_dtheta_C > 0 && myScheme.uMode == "min" || dV_dtheta_C < 0 && myScheme.uMode === "max"
                    uOpt = -uOpt
                end

                # Calculate Hamiltonian term
                result[i,j,k] = -(myScheme.obj.speed*cos_theta*dV_dx_C + myScheme.obj.speed*sin_theta* dV_dy_C + uOpt*dV_dtheta_C)


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
#end



    # Variables to calculate dissipation
    dx_dt = 0.0
    dy_dt = 0.0
    dθ_dt = 0.0
    #
    max_alpha1 = -1e9
    max_alpha2 = -1e9
    max_alpha3 = -1e9

    # Calculate dissipation
    for k = 1:myScheme.MyGrid.pts_each_dim[1] #theta index
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for i = 1:myScheme.MyGrid.pts_each_dim[3] # x index
                # Calculate alphas
                dx_dt = abs(myScheme.obj.speed * cos(thetas[k]))
                dy_dt = abs(myScheme.obj.speed * sin(thetas[k]))
                dθ_dt = abs(myScheme.obj.wMax)

                if(max_alpha1 < dx_dt)  max_alpha1 = dx_dt end
                if(max_alpha2 < dy_dt)  max_alpha2 = dy_dt end
                if(max_alpha3 < dθ_dt)  max_alpha3 = dθ_dt end

                # Calculate dissipation
                diss = 0.5*(derivDiff1[i,j,k]*dx_dt + derivDiff2[i,j,k]*dy_dt + derivDiff3[i,j,k]*dθ_dt)

                # Minus dissipipation from result
                result[i,j,k] = -(result[i,j,k] - diss)
            end
        end
    end


    # Calculate next time step
    factorCFL = 0.8

    stepBoundInv = max_alpha1/myScheme.MyGrid.dx[1]  + max_alpha2/myScheme.MyGrid.dx[2] + max_alpha3/myScheme.MyGrid.dx[3]
    stepBound = 1/stepBoundInv

    deltaT = min(factorCFL*stepBound, tspan[2] -  tspan[1])

    # Return next time point
    t = tspan[1]
    t = t + deltaT

    # Integrate V_function
    for k = 1:myScheme.MyGrid.pts_each_dim[1] # theta index
        for j = 1:myScheme.MyGrid.pts_each_dim[2] # y
            for i = 1:myScheme.MyGrid.pts_each_dim[3] # x index
                V_function[i,j,k] = V_function[i,j,k] + result[i,j,k] * deltaT
            end
        end
    end

    return  t, V_function
end
