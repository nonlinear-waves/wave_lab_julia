using DifferentialEquations

include("../bin_main/double_F.jl")

function profile_flux(p, s, s_old = nothing) 

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

    # End state tolerance
    if isnothing(s.tol)
        s_tol = 1e-4
    else
        s_tol = s.tol
    end

    # Maximum value of R allowed
    if isnothing(s.R_max)
        s_R_max = 10000
    else
        s_R_max = s.R_max
    end

    # Maximum value of L allowed
    if isnothing(s.L_max)
        s_L_max = 10000
    else
        s_L_max = s.L_max
    end

    # Numerical infinity
    s_I = 1
    # Profile solved on right half domain
    s_side = 1
    # Array for right hand side 
    s_rarray = 1:s.n
    # Array for left hand side
    s_larray = s.n + 1:(2 * s.n)

    #Bvp solver projections
    AM = s.Flinear(s.UL, p)
    proj1, _ = projection1(AM, -1, 0)
    s_LM = svd(transpose(proj1)).U

    if size(s_LM,1) == 1 && size(s_LM, 2) == 1
        n_LM = 0
    else
        n_LM = size(s_LM, 2)
    end

    AP = s.Flinear(s.UR, p)
    proj2, _ = projection1(AP, 1, 0)
    s_LP = svd(transpose(proj2)).U

    if size(s_LP,1) == 1 && size(s_LP, 2) == 1
        n_LP = 0
    else
        n_LP = size(s_LP, 2)
    end

    s_n_phs = s.n - n_LM - n_LP

    if s_n_phs < 1
        println("Eigenvalues at negative infinity: ")
        println(eigvals(AM))
        println("Eigenvalues at positive infinity")
        println(eigvals(AP))
        error("profile_flux.jl does not solve undercompressive profiles")

    end

    # BVP tolerances
    if isnothing(s.bvp_options)
        s_bvp_options = Dict(:reltol => 1e-6, :abstol => 1e-8)
    end


    new_s = ProfileSolution(s.F, s.Flinear, s.n, s.order, s.phase, s.UL, s.UR, s.stats, s_tol, s_R_max, s_L_max, s_I, nothing, nothing, s_side, s_rarray, s_larray, s_LM, s_LP, s_n_phs, s_bvp_options, nothing, nothing)

    # Positive numerical infinity

    # Solve the profile initially

    p, s = profile(p, new_s, s_old)

    # Take out extra mesh points

    # stride = how many points to take out of solution to
    # minimize points in final solution
    s_stride = 3

    return nothing

end



function profile(p, s, s_old)

    # Provided initial guess

    if !isnothing(s_old)
        if !isnothing(s_old.solver)
            if (s_old.solver != "bvp4c" || s_old.solver != "bvp5c")
                solver_type = "bvp"
                s_stride = s.stride
            else 
                solver_type = "ode"
                s_stride = 3
            end

            if solver_tpe == "ode"
                pre_guess1(x) = ode_to_bvp_guess(x, s_old, s)
                pre_guess = pre_guess1
            elseif solver_type == "bvp"
                pre_guess2(x) = continuation_guess(x, s_old, s)
                pre_guess = pre_guess2
            else
                error("Undefined solver type")
            end
        
        else 
            pre_guess3(x) = continuation_guess(x, s_old, s)
            pre_guess = pre_guess3
        end

        stride = s_old.stride
        count = 1
        for j = 1:length(s_old.sol.t)
            if stride % (j-1) == 0
                xdom[count] = s_old.sol.t[j]
                count += 1
            end
        end
        
        if (length(s_old.sol.u) - 1) % stride != 0
            x_dom[count] = s_old.sol.u[end]
        end

        s_I = s_old.I
        s_L = s_old.L
        s_R = s_old.r

    else
        s_stride = s.stride
        s_I = 1

        if isnothing(s.R)
            s_R = 5
        else
            s_R = s.R
        end

        s_L = -s_R

        pre_guess4(x) = guess(x,s)
        pre_guess = pre_guess4
        x_dom = LinRange(0,1,30)

    end

    s = ProfileSolution(s.F, s.Flinear, s.n, s.order, s.phase, s.UL, s.UR, s.stats, s.tol, s.R_max, s.L_max, s_I, s_R, s_L, s.side, s.rarray, s.larray, s.LM, s.LP, s.n_phs, s.bvp_options, s_stride, s.sol)


    # Convergence to endstates tolerance

    err = s.tol + 1
    while err > s.tol
        pre_bc(y, s, x) = bc(y, s, x)
        pre_ode(y, p_bvp, x) = double_F(y, p_bvp, x)

        p_bvp = BVP_params(s, p)



        bvp = BVProblem(pre_ode, pre_bc, pre_guess(x_dom), (x_dom[1], x_dom[end]), p_bvp)
        
        s_sol = solve(bvp, MIRK5(); s.bvp_options..., dt = 0.01)


        err1 = max(abs(s_sol.u[s.rarray, end] - s.UR))
        err2 = max(abs(s.sol.u[s.larray, end] - s.UL))
        err = maximum(err1, err2)

        if s.stats == "on"
            println("Profile Boundary error", err)
        end


        if err > s.tol
            s_R = 1.1 * s.R
            s_L = -s.R 
        end

        if err2 > s.tol
            s_L = 1.1 * s.L
            s_R = -s.L
        end

        if abs(s_L) > s.L_max
            error("Could not meet specified tolerance in profile solver without exceeding the maximum allowed value of negative infinity")
        end

        if abs(s_R) > s.R_max
            error("Could not meet specified tolerance in profile solver without exceeding the maximum allowed value of positive infinity")
        end

        if err > s.tol
            pre_guess(x) = continuation_guess(x, s_old, s)
            x_dom = s_old.sol.u
        end

        # TODO:: Figure out the correct place to update s

    end
        
    return p, s

end



function guess(x, s) 
    # Guess using tanh solution
    a = 0.5 * (s.UL + s.UR)
    c = 0.5 * (s.UL - s.UR)
    out = [a .- c * tanh.((s.R / s.I) * x) a .- c*tanh.((s.L/s.I) * x)]
    out = [out[i,:] for i in axes(out,1)]
    
    return out

end


function bc(u, p, t)
    # Boundary conditions. We split the problem in half and reflect onto the right side

    s = p.s

    #println(isa(s.order))

    if isa(s.order, Int64)
        index = s.order
    else
        index = s.order[1:s.n_phs]
    end

    out = [u[1][s.rarray] .- u[1][s.larray]
           transpose(s.LM) * u[end][s.larray] .- s.UL
           transpose(s.LP) * u[end][s.rarray] .- s.UR
           u[1][index] - s.phase[index]]


    return out

end