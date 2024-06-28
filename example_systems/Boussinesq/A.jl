function A(x, lambda, s, p)
    gamma = 0.5 * sqrt(1 - p.S^2)
    u = 1.5 * (1 - p.S^2) * sech(gamma * x) .^ 2
    ux = -2 * gamma * u .* tanh(gamma * x)
    uxx = 2 * gamma^2 * u .* (2 - 3 * sech(gamma * x) .^ 2)

    a41 = -lambda^2 - 2 * uxx
    a42 = 2 * lambda * p.S - 4 * ux
    a43 = (1 - p.S^2) - 2 * u[1,1]

    out = [0 1 0 0; 0 0 1 0; 0 0 0 1; a41 a42 a43 0]

    return out

end