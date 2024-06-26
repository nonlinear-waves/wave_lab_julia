using Test

include("../wavelab/projection.jl")
include("../example_systems/Burgers/A.jl")
include("../wavelab/structs.jl")


function test_burgers1()

    x = -20
    lambda = 1
    s = 1
    ul = 1
    ur = 0
    p = Parameter(ul, ur, "off")
    posneg = 1
    eps = 0

    actual, _ = projection2(A(x, lambda, s, p), posneg, eps)

    a = 0.5 * (ul - ur)

    mu_plus = (-a - sqrt(a^2 + 4 * lambda)) / 2
    mu_minus = (a + sqrt(a^2 + 4 * lambda)) / 2

    expected = -1 / (mu_minus - mu_plus) * [mu_plus -1
                                           mu_plus*mu_minus  -mu_minus]


    @test isapprox(actual, expected)
end


test_burgers1()