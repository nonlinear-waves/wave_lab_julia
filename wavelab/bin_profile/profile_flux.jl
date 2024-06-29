function profile_flux(p, s, s_old) 

    # Solves the profile for the flux formulation. 
    # If s.s_old exists, then this
    # is used for the initial guess. Otherwise, a tanh
    # solution is used for an initial guess. Left and right 
    # numerical infinity are expanded as needed to assure
    # the end state error is within s.tol, though s.L = s.R
    # in this program, which may not be necessary. Uneeded 
    # mesh points of the solution are removed to speed up
    # interpolation and continuation routines.
    #
    # The user must include in the input structures the following:
    #
    # s.phase - a vector of phase conditions
    # s.order - a vector indicating the order in which the phase 
    #               conditions should be applied. 
    # s.F - a function handle to the profile, e.g. s.F = @F,  F(x,y,s,p) = ...
    # s.UL, UR - end states of the n-r equations that need to be solved in
    #                   flux formulation. UL at - infinity, UR at + infinity
    # s.n = n-r in flux formulation (number of profile equations to integrate)
    #
    # Optional input:
    # s.tol - profile endstate maxim absolute error, (defaults to 1e-4)
    # s.R_max - maximum allowed interval length on right (defaults to 1000)
    # s.L_max - maximum allowed interval length on left (defaults to 1000)



    return nothing

end