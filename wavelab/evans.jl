using LinearAlgebra

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

    OmegaL0 = orth



end



