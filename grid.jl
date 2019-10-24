function add(x, y)
    return x+y
end

mutable struct grid #Must initialize first 4 elements
    grid_min::Array{Float64}
    grid_max::Array{Float64}
    N::Array{Float64}
    pDims::Int
    points_every_dimVec::Array{Float64}
    boundary_f::Array{Function}
    grid(grid_min, grid_max, N,p_dims) = new(grid_min, grid_max, N, p_dims)
end

function processGrid(grid)


end


b = add(2,3)
c =grid([1,2,3], [1,2,3],[1,2,3],3)
c.boundary_f = [add]
#d = c.boundary_f[1](2,3)
