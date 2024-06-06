using LinearAlgebra

include("manifold_polar.jl")

function evans(yl, yr, lambda, s, p, m, e)

    # Returns the evans function output at a given point.
    #
    # Input "yl" and "yr" are respectively the initializing values on the left
    # and right for the desired manifolds, "lambda" is the value in the complex
    # plane where the Evans function is evaluated, and s,p,m,e are structures
    # explained in the STABLAB documentation.


    fun = getfield(Main, Symbol(e.evans))
    
    return fun(yl, yr, lambda, s, p, m, e)

end


function reg_reg_polar(WL, WR, lambda, s, p, m, e)

    #Solve for the basis on left

    #TODO::Might want to check that this does the same thing as MATLAB code. The code uses an "econ" parameter
    OmegaL0 = svd(WL).U
    println(lambda)
    alphaL = OmegaL0' * WL
    muL = tr(OmegaL0' * e.LA(e.Li[1], lambda, s, p) * OmegaL0)
    omegal, gammal = manifold_polar(e.Li, OmegaL0, lambda, e.LA, s, p, m, e.kl, muL)

    # Use projection2 for constructing orthonormal basis

    


    return nothing
end



