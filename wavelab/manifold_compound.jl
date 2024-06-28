using LinearAlgebra
using DifferentialEquations


function manifold_compound(x, z, lambda, s, p, m, A, k, pmMU)
    # out = manifold_compound(x,z,lambda,s,p,m,A,k,pmMU)
    #
    # Returns the vector representing the manifold evaluated at x(2).
    #
    # Input "x" is the interval the manifold is computed on, "z" is the
    # initializing vector for the ode solver, "lambda" is the point on the
    # complex plane where the Evans function is computed, s,p,m are structures
    # explained in the STABLAB documentation, "A" is the function handle to the
    # desired Evans matrix, "k" is the dimension of the manifold sought, and
    # "pmMU" is 1 or -1 depending on if respectively the growth or decay 
    # manifold is sought.

    mat = A(x[1], lambda, s, p)

    D = eigvals(mat)
    R = eigvecs(mat)

    e, mat = findmax(real(pmMU * D))

    MU = D[mat]


    ode_params = ODE_params(lambda, s, p, A, m.n, k, MU)
    prob = ODEProblem(capa, x, z, ode_params)

    # According to Julia documentation, DP5 is the julia equivalent od ode45. Might want to try out Tsit5.
    # TODO:: Consider providing functionality for the user to choose the solver type
    sol = solve(prob, DP5(); m.options...)

    out = sol[end, :]

    return out

end