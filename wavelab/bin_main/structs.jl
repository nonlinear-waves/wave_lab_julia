# Numerical Infinity - typically denoted as s

struct Infinity
    I
    R
    L
    A
    Ak
end

struct ProfileSolution
    F
    Flinear
    n
    order
    phase
    UL
    UR
    stats
    tol
    R_max
    L_max
    I
    R
    L
    side
    rarray
    larray
    LM
    LP
    n_phs
    bvp_options
    stride
    sol
end


# TODO:: Figure out a better name for this struct so that it makes more sense
struct M
    n
    damping
    method
    options
    ode_fun
end


# TODO:: Figure out a better name for this struct so that it makes more sense
struct C
    LA
    RA
    stats
    refine
    tol
    ksteps
    lambda_steps
    basisL
    basisR
    evans
    epsl
    epsr
    Lproj
    Rproj
    L
    R
    check
    debug
end


# TODO:: Figure out a better name for this struct so that it makes more sense
struct E
    evans
    LA
    RA
    kl
    kr
    NL
    NR
    Li
    Ri
end

struct ODE_params
    lambda
    A
    s
    p
    n
    k
    mu
end