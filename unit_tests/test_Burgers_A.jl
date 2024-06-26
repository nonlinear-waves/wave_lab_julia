using Test

include("../example_systems/Burgers/A.jl")
include("../wavelab/structs.jl")



function test1()
    x = 1
    lambda = 1
    s = 1
    ul = 1
    ur = 0
    p = Parameter(ul, ur, "off")


    expected = [0 1
                0.8825 -0.1225]

    actual = A(x, lambda, s, p)


    @test isapprox(expected, actual, rtol = 1e-4, atol = 1e-4)

end


function test2()
    x = 2
    lambda = -1
    s = 1
    ul = 1
    ur = 0
    p = Parameter(ul, ur, "off")


    expected = [0 1
                -1.0983   -0.2311]

    actual = A(x, lambda, s, p)

    @test isapprox(expected, actual, rtol = 1e-4, atol = 1e-4)
    
end


test1()
test2()