using LinearAlgebra

function projection1(matrix,posneg,eps)
    # [P,Q] = projection1(matrix,posneg,eps)
    #
    # Returns a projector P
    #
    # Input "matrix" is the matrix from which the eigenprojection comes,
    # "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
    # sought respectively. The input eps gives a bound on how small the eigenvalues sought
    # can be, which is desirable when a zero mode should be avoided.
    
    D = eigvals(matrix)
    R = eigvecs(matrix)


    if size(R,1) == 1 && size(R,2) == 1
        R = R[0]
    end
    
    # R[:,2] = -R[:,2]


    L = inv(R)

    P = zeros(size(R))
    
    
    if posneg == 1
        index = findall(x -> x >eps, real(D))
    elseif posneg == -1
        index = findall(x -> x <eps, real(D))
    elseif posneg == 0
        index = findall(x -> x<eps, abs(real(D)))
    end
    
    for j=index
        P = P + R[:,j] .* transpose(L[j,:])
    end
    
    Q = P*R[:,index]

    return P, Q

end

function projection2(matrix, posneg, eps)

    # Returns a projector P and spanning set Q1 of the invariant subspace
    # associated with the given matrix and specified subspace.
    #
    # Input "matrix" is the matrix from which the eigenprojection comes,
    # "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
    # sought. The input eps gives a bound on how small the eigenvalues sought
    # can be, which is desirable when a zero mode should be avoided.

    # Uses Schur decomposition to get a basis for the generalized eigenspace

    #TODO: Schur decomposition does not return the same decomposition as matlab. Verify that this does not cause problems
    F = schur(matrix)

    U = F.vectors
    T = F.Schur

    # TODO:: I am slightly concerned these eigs will not be in same order as the MATLAB code. Verify that this is not the case.
    eigs = F.values

    #println(eigs)

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
        P = P + R[:, i] .* transpose(L[i, :])
    end

    return P, Q1
end


# analytic_projection([1 -2.5; 3 -0.32], 1, 0)

# function [P,Q] = projection1(matrix,posneg,eps)
#     # [P,Q] = projection1(matrix,posneg,eps)
#     #
#     # Returns a projector P
#     #
#     # Input "matrix" is the matrix from which the eigenprojection comes,
#     # "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
#     # sought respectively. The input eps gives a bound on how small the eigenvalues sought
#     # can be, which is desirable when a zero mode should be avoided.
    
#     [R,D] = eig(matrix); 
#     L = inv(R);
#     P = zeros(size(R));
    
    
#     if posneg == 1
#         index = find(real(diag(D))>eps).';
#     elseif posneg == -1
#         index = find(real(diag(D))<eps).';
#     elseif posneg == 0
#         index = find(abs(real(diag(D)))<eps).';
#     end
    
#     for j=index
#         P = P + R(:,j)*L(j,:);
#     end
    
#     Q = P*R(:,index);