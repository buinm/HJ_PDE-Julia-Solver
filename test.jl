

#g = makeGrid([-5.0,-5.0,-pi], [5,5,pi],[41,41,41],3)
#my_car = makeDubinsCar([0.0,0.0,0.0],1.0, 1.0)
#myScheme = sf(MyGrid = g, obj=my_car)


a = zeros(40,40,40)
b = zeros(40,40,40)
b = b.*10

c = zeros(40,40,40)
d = zeros(40,40,40)


@time begin
 for i = 1:40
    for j = 1:40
         for k = 1:40
            d[i,j,k] = a[i,j,k] + b[i,j,k]
        end
    end
end
end


@time begin
c = a + b
end
