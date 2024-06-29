using DifferentialEquations

function manifold_polar(x, y, lambda, A, s, p, m, k, mu)
    # Returns "Omega", the orthogonal basis for the manifold evaluated at x(2)
    # and "gamma" the radial equation evaluated at x(2).
    #
    # Input "x" is the interval on which the manifold is solved, "y" is the
    # initializing vector, "lambda" is the point in the complex plane where the
    # Evans function is evaluated, "A" is a function handle to the Evans
    # matrix, s, p,and m are structures explained in the STABLAB documentation, 
    # and k is the dimension of the manifold sought.

    # Solve the ODE 

    #unused, Y = m.ode_fun(m.method, [x[1] x[2]])

    ode_params = ODE_params(lambda, A, s, p, m.n, k, mu)


    prob = ODEProblem(m.method, [reshape(y, m.n*k, 1); 0], (x[1], x[2]), ode_params)

    #TODO:: It doesn't look like drury works with automatic differentiation. The issue has to do with Complex types? Look into modifying drury so it is compatible
    sol = solve(prob, m.ode_fun(autodiff=false); m.options...)

    #TODO:: This indexing is a little funky and might cause problems for higher dimensional problems
    Omega = reshape(sol[1:m.n * k, :, end], m.n, k);

    gamma = exp(sol[m.n*k + 1, :, end][1])

    return Omega, gamma
end