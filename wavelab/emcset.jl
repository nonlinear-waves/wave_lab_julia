using DifferentialEquations

include("Aadj.jl")
include("Akadj.jl")
include("analytic_basis.jl")
include("evans.jl")
include("projection.jl")
include("drury.jl")
include("structs.jl")


# TODO:: Fix the default arguments for func and compound_func to reflect MATLAB behavior

function emcset(s, shock_type, eLR, Evan_type = "default", func = nothing, compound_func = nothing)
    
    #     function [e,m,c] = emcset(shock_type,eL,eR,Evan_type)
    #
    # Sets the values of the STABLAB structures e, m, and c to 
    # default values. Takes as input a string, shock_type, which is either
    # "front" or "periodic". The input eL and eR are respectively the 
    # dimension of the left and right eigenspaces of the Evans matrix.
    # The input Evan_type is an optional string. If not specified, Evan_type
    # will be assigned the most advantageous polar coordinate method.
    # Evan_type has the following options when shock_type = 'front':
    # 
    # reg_reg_polar
    # reg_adj_polar
    # adj_reg_polar
    # reg_adj_compound
    # adj_reg_compound
    # 
    # when shock_type = 'periodic', the choices are:
    # 
    # regular_periodic
    # balanced_periodic
    # balanced_polar_scaled_periodic
    # balanced_polar_periodic
    # balanced_scaled_periodic


    eL = eLR[1];
    eR = eLR[2];

    # TODO:: Figure out how to make the func and compound_func default to A and Ak where those are files that are in the current directory or path like in MATLAB

    if cmp(shock_type, "front") == 0
        e, m, c = initialize_front(s, eL, eR, Evan_type, func, compound_func)
        new_s = Infinity(s.I, s.R, s.L, func, compound_func)
    
    else
        error("User must specify which type of traveling wave is being studied")

    
    end

    # TODO:: Modify return function when function is finished
    return new_s, e, m, c
end


function initialize_front(s, kL, kR, Evan_type, func, compound_func)

    m_n = kL + kR

    if cmp(Evan_type, "default") == 0

        if kL > m_n / 2
            e_evans = "adj_reg_polar"

        elseif kL < m_n / 2
            e_evans = "reg_adj_polar"

        else
            e_evans = "reg_reg_polar"

        end

    else
        e_evans = Evan_type
    end

    println(e_evans)
    if cmp(e_evans, "reg_adj_polar") == 0
        c_LA = func
        e_LA = c_LA
        c_RA = Aadj
        e_RA = c_RA
        e_kl = kL
        e_kr = n - kR

        # TODO:: Verify that this won't cause problems elsewhere (NOT in MATLAB code)
        e_NL = 0
        e_NR = 0

    elseif cmp(e_evans, "reg_reg_polar") == 0
        println(func)
        c_LA = func
        e_LA = c_LA
        c_RA = func
        e_RA = c_RA
        e_kl = kL
        e_kr = kR

        # TODO:: Verify that this won't cause problems elsewhere (NOT in MATLAB code)
        e_NL = 0
        e_NR = 0

    elseif cmp(e_evans, "adj_reg_polar") == 0
        c_LA = Aadj
        e_LA = Aadj
        c_RA = func
        e_RA = func
        e_kl = n - kL
        e_kr = kR

        # TODO:: Verify that this won't cause problems elsewhere (NOT in MATLAB code)
        e_NL = 0
        e_NR = 0

    elseif cmp(e_evans, "reg_adj_compound") == 0
        c_LA = func
        e_LA = compound_func
        c_RA = Aadj
        e_RA = Akadj
        e_kl = kL
        e_kr = kR

        # TODO:: Verify that this won't cause problems elsewhere (NOT in MATLAB code)
        e_NL = 0
        e_NR = 0

    elseif cmp(e_evans, "adj_reg_compound") == 0
        c_LA = Aadj
        e_LA = Akadj
        c_RA = func
        e_RA = compound_func
        e_kl = kL
        e_kr = kR

        # TODO:: Verify that this won't cause problems elsewhere (NOT in MATLAB code)
        e_NL = 0
        e_NR = 0

    elseif cmp(e_evans, "reg_reg_cheby") == 0
        c_LA = func
        e_LA = c_LA
        c_RA = func
        e_RA = c_RA
        e_kl = kL
        e_kr = kR
        e_NL = 60
        e_NR = 60

    end

    c_stats = "off"
    c_refine = "off"
    c_tol = 0.2
    c_ksteps = 2^5
    c_lambda_steps = 0
    c_basisL = analytic_basis
    c_basisR = analytic_basis
    c_evans = evans

    c_epsl = 0
    c_epsr = 0
    c_Lproj = projection2
    c_Rproj = projection2

    m_damping = 0
    m_method = drury

    # TODO:: Not sure what to do with the options parameters yet
    m_options = 0

    # TODO:: Decide which Julia solver would be best here, or if we should include functionality for the user to pick
    m_ode_fun = QNDF

    # Dependent structure variables
    e_Li = [s.L 0]
    e_Ri = [s.R 0]
    c_L = s.L
    c_R = s.R

    #Create structures
    m = M(m_n, m_damping, m_method, m_options, m_ode_fun)
    c = C(c_LA, c_RA, c_stats, c_refine, c_tol, c_ksteps, c_lambda_steps, c_basisL, c_basisR, c_evans, c_epsl, c_epsr, c_Lproj, c_Rproj, c_L, c_R)
    e = E(e_evans, e_LA, e_kl, e_kr, e_NL, e_NR)


    # TODO:: Modify return function when function is finished
    return e, m, c

end

    