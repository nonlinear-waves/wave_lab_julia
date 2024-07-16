function double_F(y, p, x, args...) 
    # Returns the split domain for the ode given in the function F.
    #
    # Input "x" and "y" are provided by the ode solver.Note that s.rarray
    # should be [1,2,...,k] and s.larray should be [k+1,k+2,...,2k]. See
    # STABLAB documentation for more inforamtion about the structure s.


    # TODO:: Not quite sure what to do about the variable number of arguments
    # Think ... could be a solution or just defining 
    
    out = [(p.s.R / p.s.I) * p.s.F((p.s.R/p.s.I) * x, y[p.s.rarray, :] , p.s, p.p, args...); (p.s.L/p.s.I) * p.s.F((p.s.L/p.s.I) * x, y[p.s.larray, :], p.s, p.p, args...)]
    return out
end