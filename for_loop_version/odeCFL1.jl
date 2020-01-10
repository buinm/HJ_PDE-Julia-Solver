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
    end

@time begin
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
end

function odeCFL1(tspan, V_function, myScheme, result, derivDiff1, derivDiff2, derivDiff3)
    # Calculate the Hamiltonian terms
    min_deriv_1 = 1e9
    min_deriv_2 = 1e9
    min_deriv_3 = 1e9
    max_deriv_1 = -1e9
    max_deriv_2 = -1e9
    max_deriv_3 = -1e9

    #result = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    #derivDiff1 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    #derivDiff2 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    #derivDiff3 = Array{Float64, 3}(undef, myScheme.MyGrid.pts_each_dim[1],myScheme.MyGrid.pts_each_dim[2],myScheme.MyGrid.pts_each_dim[3])
    dV_dx_L =  0.0
    dV_dx_R = 0.0
    dV_dy_L = 0.0
    dV_dy_R = 0.0
    dV_dtheta_L = 0.0
    dV_dtheta_R = 0.0
    uOpt = 0.0
    thetas = myScheme.MyGrid.min[3] : myScheme.MyGrid.dx[3] : myScheme.MyGrid.max[3]
    cos_theta = 0.0
    sin_theta = 0.0

    #@time begin
    for k = 1:myScheme.MyGrid.pts_each_dim[3] #theta index
            theta = thetas[k]
            cos_theta = cos(theta)
            sin_theta = sin(theta)
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for i = 1:myScheme.MyGrid.pts_each_dim[1] # x
                # Left and Right spatial derivative
                #@time begin
                @views dV_dx_L, dV_dx_R = spatial_Derivative(myScheme.MyGrid, i, V_function[:, j,k], 1)
            #end
                #@time begin
                @inbounds @views dV_dy_L, dV_dy_R = spatial_Derivative(myScheme.MyGrid, j, V_function[i, :,k], 2)
                #end
                @inbounds @views dV_dtheta_L, dV_dtheta_R = spatial_Derivative(myScheme.MyGrid, k, V_function[i, j,:], 3)

                # Average Spatial Derivative
                dV_dx_C = (dV_dx_L + dV_dx_R)/2
                dV_dy_C = (dV_dy_L + dV_dy_R)/2
                dV_dtheta_C = (dV_dtheta_L + dV_dtheta_R)/2

                # Calculate optimal control input
                uOpt = myScheme.obj.speed
                if (dV_dtheta_C > 0 && myScheme.uMode === "min") || (dV_dtheta_C < 0 && myScheme.uMode === "max")
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

    # Calculate dissipation
    dx_dt = 0.0
    dy_dt = 0.0
    dθ_dt = 0.0
    #
    max_alpha1 = -1e9
    max_alpha2 = -1e9
    max_alpha3 = -1e9

    for k = 1:myScheme.MyGrid.pts_each_dim[1] #theta index
        for j = 1:myScheme.MyGrid.pts_each_dim[2] #y
            for i = 1:myScheme.MyGrid.pts_each_dim[3]
                # Calculate alphas
                dx_dt = abs(myScheme.obj.speed * cos(thetas[k]))
                dy_dt = abs(myScheme.obj.speed * sin(thetas[k]))
                dθ_dt = abs(myScheme.obj.wMax)

                if(max_alpha1 < dx_dt)  max_alpha1 = dx_dt end
                if(max_alpha2 < dy_dt)  max_alpha2 = dy_dt end
                if(max_alpha3 < dθ_dt)  max_alpha3 = dθ_dt end

                # Calculate dissipation
                diss = 0.5*(derivDiff1[i,j,k]*dx_dt + derivDiff2[i,j,k]*dy_dt + derivDiff3[i,j,k]*dθ_dt)

                # I
                result[i,j,k] = diss -result[i,j,k]
            end
        end
    end

    # Calculate next time step
    stepBoundInv = 0
    stepBoundInv = max_alpha1/myScheme.MyGrid.dx[1]  + max_alpha2/myScheme.MyGrid.dx[2] + max_alpha3/myScheme.MyGrid.dx[3]
    delta_t = 1/stepBoundInv

    # Return next time point
    t = tspan[1]
    t = t + delta_t

    # Integrate V_function
    V_function = V_function .+ result .* delta_t

    return  t, V_function
end
