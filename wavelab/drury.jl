function drury(y, p, t)
    # Returns the ODE output for the polar method using the method of Drury
    #
    # Input "t" and "y" are provided by the ode soler, "A" is a function handle to the
    # desired Evans matrix, s and p are structures explained in the STABLAB
    # documentation, "n" is the dimension of the system and "k" is the
    # dimension of the manifold

    # W = Omega*alpha, gamma = det(alpha), rho = log(gamma)

    # Reshape the (k*n + 1) vector to be the n x k matrix Omega
    Omega = reshape(y[1:p.k*p.n,1], p.n, p.k)

    # Evaluate A(x, lamda)
    A_temp = A(t, p.lambda, p.s, p.p)

    # Compute Omega' and rho'
    return [reshape((I - Omega * Omega') * A_temp * Omega, p.n * p.k, 1); tr(Omega' * A_temp * Omega) - p.mu]

end