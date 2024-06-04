

# Numerical Infinity - typically denoted as s

struct Infinity
    I
    R
    L
    A
    Ak
end

# Parameters - typically denoted as p
struct Parameter
    ul
    ur
    integrated
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
end


# TODO:: Figure out a better name for this struct so that it makes more sense
struct E
    evans
    LA
    kl
    kr
    NL
    NR
end