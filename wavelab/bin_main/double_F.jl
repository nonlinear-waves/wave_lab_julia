function double_F(y, p, x, args...) 
    # Returns the split domain for the ode given in the function F.
    #
    # Input "x" and "y" are provided by the ode solver.Note that s.rarray
    # should be [1,2,...,k] and s.larray should be [k+1,k+2,...,2k]. See
    # STABLAB documentation for more inforamtion about the structure s.


    # TODO:: Not quite sure what to do about the variable number of arguments
    # Think ... could be a solution or just defining 
    s = p.s
    p = p.p

    out = [(s.R / s.I) * s.F((s.R/s.I) * x, y[s.rarray, :] , s, p, args...); (s.L/s.I) * s.F((s.L/s.I) * x, y[s.larray, :], s, p, args...)]

    return out
end