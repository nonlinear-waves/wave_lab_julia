
function analytic_basis(projection, x, preimage, s, p, A, posneg, eps, Q = nothing, p_old_in = nothing)
    # Returns an analytic basis at specified infinity using the method of Kato.
    # 
    # Input "projection" is a function handle to the projection function to be
    # used, "x" is the numerical value of infinity, "preimage" is the contour
    # on which the Evans function is computed, "s" and "p" are structures
    # explained in the STABLAB documentation, "A" is a function handle to the
    # Evans matrix, "posneg" is 1 or -1 determining which space the
    # projection function should return, and "eps" is the tolerance in the
    # projection function. If input Q and p_old_in are specified, then the
    # first value of the analtyic basis is set to Q with projection p_old_in.
    # This allows for continuing or filling in a previous analtyic basis
    # computation. When these two inputs are not specified, they are chosen
    # automtically.

    #TODO:: There may be problems with this code if projection is not actually one dimensional. Verify if this will ever be the case
    iterations = size(preimage, 1)


    # println(projection)
    # println(x)
    # println(s)
    # println(p)
    # println(A)
    # println(posneg)
    # println(eps)

    p_old, Q1 = projection(A(x, preimage[1], s, p), posneg, eps)
    n,k = size(Q1)
    out = zeros(ComplexF64, n, k, iterations)

    projects = zeros(ComplexF64, size(p_old, 1), size(p_old, 2), iterations)

    if isnothing(Q)
        out[:,:,1] = Q1
        projects[:,:,1] = p_old

    else
        out[:,:,1] = Q
        projects[:,:,1] = p_old_in
        p_old = p_old_in

    end

    for j=2:iterations
        proj, _ = projection(A(x, preimage[j], s, p), posneg, eps)

        # a = (p.ul - p.ur) / 2
        # mu_plus = (-a - sqrt(a^2 + 4 * preimage[j])) / 2
        # mu_minus = (a + sqrt(a^2 + 4 * preimage[j])) / 2

        # println("mu_plus:", (-a - sqrt(a^2 + 4 * 1)) / 2)
        # println("mu_minus:", (a + sqrt(a^2 + 4 * 1)) / 2)

        # proj = 1 / (mu_minus - mu_plus) * [mu_minus 1
        #                                    mu_minus * mu_plus   -mu_plus]
        # println(proj)
        out[:,:,j] = proj * (I + 0.5 * p_old * (I - proj)) * out[:, :, j-1]
        projects[:,:,j] = proj
        p_old = proj

    end

    # println(out)

    return out, projects


end

