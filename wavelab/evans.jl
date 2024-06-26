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
    alphaL = OmegaL0' * WL
    muL = tr(OmegaL0' * e.LA(e.Li[1], lambda, s, p) * OmegaL0)
    omegal, gammal = manifold_polar(e.Li, OmegaL0, lambda, e.LA, s, p, m, e.kl, muL)
    # println(e.Li)
    # println(OmegaL0)
    # println(lambda)
    # println(e.LA)
    # println(s)
    # println(p)
    # println(m)
    # println(e.kl)
    # println(muL)
    # println("OmegaL", omegal)
    # println("gammaL", gammal)

    # Solve for the basis on the right

    OmegaR0 = svd(WR).U
    alphaR = OmegaR0' * WR
    muR = tr(OmegaR0' * e.RA(e.Ri[1], lambda, s, p) * OmegaR0)
    omegar, gammar = manifold_polar(e.Ri, OmegaR0, lambda, e.RA, s, p, m, e.kr, muR)
    # println(e.Ri)
    # println(OmegaR0)
    # println(lambda)
    # println(e.RA)
    # println(s)
    # println(p)
    # println(m)
    # println(e.kr)
    # println(muR)
    # println("omegaR", omegar)
    # println("gammar", gammar)

    #Evaluate the determinant

    out = (det(alphaL) * gammal) * (det(alphaR) * gammar) * det([omegal omegar])

    return out
    
end



