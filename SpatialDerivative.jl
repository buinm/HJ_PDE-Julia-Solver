include("./grid.jl")
include("./BasicShapes.jl")
include("./BoundaryCondition.jl")
using .Grid
using .BasicShapes
using .BoundaryCondition



function upwindFirstFirst(grid, data, dim)
    # Will add ghost function later, work on data array

    gdata = Array{Float64, grid.dim}(undef, tuple(grid.pts_each_dim...,))
    if dim != grid.pDims
        gdata = addGhostExtrapolate(data,dim,1,true)
    else
        gdata = addGhostPeriodic(data,dim,1,false)
    end

    sizeData = size(gdata)
    indices1 = Array{Any, 1}(undef, grid.dim)
    indices2 = Array{Any, 1}(undef, grid.dim)

    for i = 1:grid.dim
        if i == dim
            indices1[i] = 2:sizeData[i]
            indices2[i] = 1:sizeData[i]-1
        else
            indices1[i] = 1:sizeData[i]
            indices2[i] = indices1[i]
        end

    end

    cartesianInd1 = CartesianIndices(tuple(indices1...,))
    cartesianInd2 = CartesianIndices(tuple(indices2...,))

    dxInv = 1/grid.dx[dim]
    deriv = dxInv*(gdata[cartesianInd1] - gdata[cartesianInd2])


    # Take leftmost data for left approximation
    indices1[dim] = 1:size(deriv,dim) - 1
    indices2[dim] = 2:size(deriv,dim)

    cartesianInd1 = CartesianIndices(tuple(indices1...,))
    cartesianInd2 = CartesianIndices(tuple(indices2...,))

    derivL = deriv[cartesianInd1]
    derivR = deriv[cartesianInd2]

    return derivL, derivR
    # Take rightmost data for right approximation
end

my_grid = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
my_data = CyclinderShape(my_grid, 3.0,[0;0;0],1.0)
l,r = upwindFirstFirst(my_grid, my_data,2)
