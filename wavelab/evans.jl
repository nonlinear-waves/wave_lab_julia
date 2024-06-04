

function evans(yl, yr, lambda, s, p, m, e)

    # Returns the evans function output at a given point.
    #
    # Input "yl" and "yr" are respectively the initializing values on the left
    # and right for the desired manifolds, "lambda" is the value in the complex
    # plane where the Evans function is evaluated, and s,p,m,e are structures
    # explained in the STABLAB documentation.


    # TODO:: Original MATLAB code uses str2func here. This might cause problems in the future
    fun = e_evans
    
    return fun(yl, yr, lambda, s, p, m, e)

end



