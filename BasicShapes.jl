#module BasicShape
include("./grid.jl")
using .Grid


function CyclinderShape(grid,ignoreDims, center, radius)
    data = zeros(tuple(grid.pts_each_dim ...,))

    #ns_per_dim = [10;10;12]
    ns_each_dim = copy(grid.pts_each_dim)
    for i = 1:grid.dim
        tmp = ns_each_dim[1]
        ns_each_dim[1] = ns_each_dim[i]
        ns_each_dim[i] = tmp
        data = reshape(data,tuple(ns_each_dim ...,))
        if i != grid.pDims
            #Broadcastable data
            data = data .+ (grid.vs[i] .- center[i]).^2
        end
    end

    data = sqrt.(data) .- radius
    return reshape(data,tuple(grid.pts_each_dim ...,))
end

c = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[40,41,42],3)
my_data = CyclinderShape(c, 3.0,[0;0;0],1.0)
#println("Result is: ",my_data)
