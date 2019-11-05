module BasicShapes

include("./grid.jl")
using .Grid

export CyclinderShape

function CyclinderShape(grid,ignoreDims, center, radius)
    data = zeros(tuple(grid.pts_each_dim ...,))

    for i = 1:grid.dim
        if i != ignoreDims
            data = data .+ (grid.vs[i] .- center[i]).^2
        end
    end
    data = sqrt.(data) .- radius
    return data
end

#c = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
#my_data = CyclinderShape(c, 3.0,[0;0;0],1.0)
#println("Result is: ",my_data)
end
