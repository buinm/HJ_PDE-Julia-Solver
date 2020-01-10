include("./DubinsCar.jl")
include("./for_loop_version\\for_grid.jl")


#g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
#my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0)
#myScheme = sf(MyGrid = g, obj=my_car)

let
    #result = Array{Float64, 3}(undef, 40,40,40)
    #cos_theta = 100
    #@time begin
    #    result[3,3,3] = cos_theta
    #end

    g = makeNewGrid([-5.0,-5.0,-pi], [5,5,pi],[40,40,40],3)
    my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0, [0,0,0])
    result = 0
    thetas = g.min[3] : g.dx[3] : g.max[3]
    println(typeof(thetas))
    println(thetas)
        for j = 1:5
            @time begin
                result = cos(thetas[j])*my_car.speed
            end
    end
end
