module Grid

# Add libary using
import Base.Iterators
export grid, makeGrid#addGhostPeriodic, addGhostExtrapolate # export module to be used by others

# TODO: implement these 2 below functions (or maybe include it from another files?)


mutable struct grid #Must initialize first 4 elements
    min::Vector{Float64}
    max::Vector{Float64}
    pts_each_dim::Vector{Int}
    pDims::Int
    # This below is partial constructor
    grid(min, max, pts_each_dim,pDims) = new(min, max, pts_each_dim, pDims)
    dim::Int
    dx::Vector{Float64}
    vs::Array{Any, 1}
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
    g.max[g.pDims] = g.min[g.pDims] + (g.max[g.pDims]-g.min[g.pDims]) * (1-1/g.pts_each_dim[g.pDims])
    g.dx =  (g.max - g.min)./(g.pts_each_dim .- 1.0)

    #g.vs = Array{Vector{Float64},1}(undef,g.dim)
    g.vs = Array{Any, 1}(undef, g.dim)
    #g.xs = Array{Any, g.dim}(Tuple(g.pts_each_dim))
    # Infer positions in each dimension
    for i = 1:g.dim
        g.vs[i] = collect(range(g.min[i], step = g.dx[i], stop = g.max[i]))
    end

    for i = 1:g.dim
        broadcast_map = ones(Int,g.dim)
        broadcast_map[i] = g.pts_each_dim[i]
        g.vs[i] = reshape(g.vs[i], tuple(broadcast_map...,))
    end

end

function makeGrid(gridMin, gridMax, num_g_pts, pDims)
    g = grid(gridMin, gridMax, num_g_pts, pDims)
    processGrid(g)
    return g
end

end
