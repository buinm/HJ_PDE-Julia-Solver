module BoundaryCondition

include("./grid.jl")
using .Grid

export addGhostPeriodic, addGhostExtrapolate


function addGhostExtrapolate(dataIn, dim, width, towardZero)
    slopeMultiplier = 0
    if towardZero
        #println("TowardZero")
        slopeMultiplier = -1
    else
        #println("Not TowardZero")
        slopeMultiplier = +1
    end
    dims = ndims(dataIn)
    my_size = size(dataIn)
    sizeIn = Array{Any, 1}(undef, dims)
    for k = 1:dims
        sizeIn[k] = my_size[k]
    end
    println(sizeIn)

    indicesIn = Array{Any, 1}(undef, dims)
    indicesOut = Array{Any, 1}(undef, dims)


    for i = 1:dims
        indicesIn[i] = 1:sizeIn[i]
    end
    indicesOut = copy(indicesIn)

    # Create sized output array
    sizeOut = copy(sizeIn)
    sizeOut[dim] = sizeIn[dim] + 2*width
    dataOut = zeros(tuple(sizeOut...,))

    # Fill in middle part of dataOut with dataIn
    indicesOut[dim] = (width+1):(sizeOut[dim] -width)
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    dataOut[cartesianIndOut] = copy(dataIn)

    # Compute Bottom Slope
    indicesOut[dim] = 1
    indicesIn[dim]  = 2
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))

    slopeBot =  dataIn[cartesianIndOut]-dataIn[cartesianIndIn]

    # Compute Top Slope
    indicesOut[dim] = sizeIn[dim]
    indicesIn[dim]  = sizeIn[dim] - 1
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
    slopeTop = dataIn[cartesianIndOut]-dataIn[cartesianIndIn]

    # Adjust slope values with sign of data at array edge
    indicesIn[dim] =  1
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
    slopeBot = slopeMultiplier * abs.(slopeBot) .* sign.(dataIn[cartesianIndIn])
    indicesIn[dim] =  sizeIn[dim]
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
    slopeTop = slopeMultiplier * abs.(slopeTop) .* sign.(dataIn[cartesianIndIn])

    # Extrapolate
    for i = 1:width
        indicesOut[dim] = i
        indicesIn[dim]  = 1
        cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
        cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
        dataOut[cartesianIndOut] = dataIn[cartesianIndIn] + (width - i+1)* slopeBot

        indicesOut[dim] = sizeOut[dim] - i + 1
        indicesIn[dim]  = sizeIn[dim]
        cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
        cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
        dataOut[cartesianIndOut] = dataIn[cartesianIndIn] + (width - i+1)* slopeTop
    end

    return dataOut
end

function addGhostPeriodic(dataIn, dim, width, towardZero)
    slopeMultiplier = 0
    if towardZero
        slopeMultiplier = -1
    else
        slopeMultiplier = +1
    end
    dims = ndims(dataIn)
    my_size = size(dataIn)
    sizeIn = Array{Any, 1}(undef, dims)
    for k = 1:dims
        sizeIn[k] = my_size[k]
    end
    # These are used to index into data array
    indicesIn = Array{Any, 1}(undef, dims)
    indicesOut = Array{Any, 1}(undef, dims)

    for i = 1:dims
        indicesIn[i] = 1:sizeIn[i]
    end
    indicesOut = copy(indicesIn)
    #println("I'm here")
    # Create sized output array
    sizeOut = copy(sizeIn)
    sizeOut[dim] = sizeIn[dim] + 2*width
    dataOut = zeros(tuple(sizeOut...,))
    #println(size(dataOut))
    #println("I'm here")

    # Fill in middle part of dataOut with dataIn
    indicesOut[dim] = width+1:sizeOut[dim] -width
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    dataOut[cartesianIndOut] = dataIn
#    println(dataOut)

    #println("I'm here")

    # Fill ghost cells
    indicesIn[dim] = (sizeIn[dim]-width+1):sizeIn[dim]
    indicesOut[dim] = 1:width
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
    dataOut[cartesianIndOut] = dataIn[cartesianIndIn]
    #println(dataOut)


    indicesIn[dim] = 1:width
    indicesOut[dim] = (sizeOut[dim]-width+1):sizeOut[dim]
    cartesianIndOut = CartesianIndices(tuple(indicesOut...,))
    cartesianIndIn  = CartesianIndices(tuple(indicesIn...,))
    dataOut[cartesianIndOut] = dataIn[cartesianIndIn]

    return dataOut


end

end
