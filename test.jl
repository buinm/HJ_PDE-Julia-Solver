
function testArray(A)
    for i = 1:4
        A[i] = 5
    end
end

B = zeros(10)
testArray(B)
