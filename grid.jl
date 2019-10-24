#module Grid

# Add libary using
#using LinearAlgebra
#export grid, processGrid, addGhostPeriodic, addGhostExtrapolate # export module to be used by others

# TODO: implement these 2 below functions (or maybe include it from another files?)
function addGhostPeriodic(dataIn, dim, width, towardZero)
    slopeMultiplier = 0
    if towardZero
        slopeMultiplier = -1
    else
        slopeMultiplier = +1
    end
    dims = ndims(dataIn)
    sizeIn = size(dataIn)
end

function addGhostExtrapolate()

end


mutable struct grid #Must initialize first 4 elements
    min::Vector{Float64}
    max::Vector{Float64}
    pts_each_dim::Vector{Float64}
    pDims::Int
    # This below is partial constructor
    grid(min, max, pts_each_dim,pDims) = new(min, max, pts_each_dim, pDims)
    dim::Int
    dx::Vector{Float64}
    vs::Array{Any, 1}
    xs::Vector{Float64}
end

# NOTE: g is a mutable struct, hence changes made here are valid
function processGrid(g) # This should turn g into a complete grid structure
    ######################
    # Output grid now should include the following fields:
    # -  vs
    # -  dx
    # -  xs
    ######################

    # Infer system dimensions
    g.dim = length(g.min)
    g.dx =  (g.max - g.min)./(g.pts_each_dim .- 1.0)
    #g.vs = Array{Vector{Float64},1}(undef,g.dim)
    g.vs = Array{Any, 1}(undef, g.dim)
    # Infer positions in each dimension
    for i = 1:g.dim
        g.vs[i] = collect(range(g.min[i], step = g.dx[i], stop = g.max[i]))
    end

    g.max[g.pDims] = g.min[g.pDims] + (g.max[g.pDims]-g.min[g.pDims]) * (1-1/g.pts_each_dim[g.pDims])
end



#b = add(2,3)
c =grid([1,2,3], [4,5,6],[10,10,12],3)
#c.bdry_f = [add]
processGrid(c)

#end
