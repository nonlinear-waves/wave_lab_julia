using DifferentialEquations

include("../../wavelab/structs.jl")
include("../../wavelab/emcset.jl")
include("../../wavelab/contour.jl")
include("../../wavelab/winding_number.jl")


include("A.jl")
include("Ak.jl")


# Parameters - typically denoted as p
struct Parameter
    B
    gamma
    mu0
    sigma
    vp
    mu
    eta
    a   
end


function MHD() 

    # Parameters

    p_B = 2
    p_gamma = 5/3
    p_mu0 = 1
    p_sigma = 1
    p_vp = 0.0001
    p_mu = 1
    p_eta = -2 * p_mu /3

    # Dependent Parameters

    p_a = p_vp ^ p_gamma * ((1 - p_vp) / (1 - p_vp ^ p_gamma))

    p = Parameter(p_B, p_gamma, p_mu0, p_sigma, p_vp, p_mu, p_eta, p_a)


    # Profile Solution


    s_F(x, y, s, p) = (2 * p.mu + p.eta) ^ (-1) * y * (y - 1 + p.a * (y ^ (-p.gamma) - 1))
    s_Flinear(y, p) = (2 * p.mu + p.eta) ^ (-1) * (y - 1 + p.a * (y ^ (-p.gamma) - 1)) + (2 * p.mu + p.eta) ^ (-1) * y * (1 - p.a * p.gamma * y ^ (-p.gamma - 1))
    # number of profile quations to integrate
    s_n = 1
    s_order = 1 
    s_phase = 0.5 * (1 + p_vp)
    # end states
    s_UL = 1
    s_UR = p_vp

    s_stats = "on"
    #Tolerance at end states
    s_tol = 1e-6

    s = ProfileSolution(s_F, s_Flinear, s_n, s_order, s_phase, s_UL, s_UR, s_stats, s_tol)

    p, s = profile_flux(p,s)


    # Structure Variables
    s, e, m, c = emcset(s, "front", [2,2], "default", A, Ak)  #default for MHD1 is reg_reg_polar
    # s, e, m, c = emcset(s, "front", [2,2], "reg_adj_polar", A)
    # s, e, m, c = emcset(s, "front", [2,2], "adj_reg_polar", A)
    # s, e, m, c = emcset(s, "front", [2,2], "reg_reg_polar", A)
    # s, e, m, c = emcset(s, "front", [2,2], "reg_adj_compound", A, Ak)
    # s, e, m, c = emcset(s, "front", [2,2], "adj_reg_compound", A , Ak)


    m_ode_fun = QNDF


    # Preimage contour

    R = 1
    circpnts = 30
    imagpnts = 30
    spread = 4
    zerodist = 1e-4


    preimage = semicirc(circpnts, imagpnts, c.ksteps, R, spread, zerodist)


    # Compute Evans function

    halfw = contour(c, s, p, m, e, preimage)
    w = [halfw; reverse(conj(halfw))]

    # Process and display data 

    wnd = 


end


MHD()