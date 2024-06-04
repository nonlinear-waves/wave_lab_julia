using LinearAlgebra

function projection2(matrix, posneg, eps)

    # Returns a projector P and spanning set Q1 of the invariant subspace
    # associated with the given matrix and specified subspace.
    #
    # Input "matrix" is the matrix from which the eigenprojection comes,
    # "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
    # sought. The input eps gives a bound on how small the eigenvalues sought
    # can be, which is desirable when a zero mode should be avoided.

    # Uses Schur decomposition to get a basis for the generalized eigenspace

    F = Schur{Complex}(schur(matrix))

    U = F.vectors
    T = F.Schur

    println(U)

    # TODO:: I am slightly concerned these eigs will not be in same order as the MATLAB code. Verify that this is not the case.
    eigs = F.values

    signed_eigs = posneg * real(eigs)
    k = length(signed_eigs[signed_eigs .> eps])
    FS = ordschur(F, signed_eigs .> eps)
    US = FS.vectors

    Q1 = US[:, 1:k]

    F_neg = Schur{Complex}(schur(-matrix))
    
    U_neg = F_neg.vectors
    T_neg = F_neg.Schur
    eigs_neg = F_neg.values

    signed_eigs_neg = posneg * real(eigs_neg)
    k_neg = length(signed_eigs_neg[signed_eigs_neg .> -eps])
    FS_neg = ordschur(F_neg, signed_eigs_neg .> -eps)
    US_neg = FS_neg.vectors

    Q2 = US_neg[:, 1:k_neg]

    R = [Q1 Q2]

    L = inv(R)

    P = zeros(size(matrix))

    for i = 1:size(Q1,2)
        P = P + R[:, i] * L[i, :]
    end

    println(P * U)


    return P, Q1
end