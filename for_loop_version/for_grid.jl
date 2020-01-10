#module Grid

# Add libary using
import Base.Iterators
#export grid, makeGrid#addGhostPeriodic, addGhostExtrapolate # export module to be used by others

# TODO: implement these 2 below functions (or maybe include it from another files?)


mutable struct new_grid #Must initialize first 4 elements
    min::Vector{Float64}
    max::Vector{Float64}
    pts_each_dim::Vector{Int}
    pDims::Int
    # This below is partial constructor
    new_grid(min, max, pts_each_dim,pDims) = new(min, max, pts_each_dim, pDims)
    dim::Int
    dx::Vector{Float64}
    vs::Array{Any, 1}
end

# NOTE: g is a mutable struct, hence changes made here are valid
function processNewGrid(g) # This should turn g into a complete grid structure
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

    g.vs = Array{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, 1}(undef, g.dim)

    # Infer positions in each dimension
    for i = 1:g.dim
        g.vs[i] = g.min[i] : g.dx[i] : g.max[i]
    end

end

function makeNewGrid(gridMin, gridMax, num_g_pts, pDims)
    g = new_grid(gridMin, gridMax, num_g_pts, pDims)
    processNewGrid(g)
    return g
end

#end
