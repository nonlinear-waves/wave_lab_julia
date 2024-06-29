using Plots


struct Parameter
    S
end

include("A.jl")
include("Ak.jl")

include("../../wavelab/structs.jl")
include("../../wavelab/emcset.jl")
include("../../wavelab/contour.jl")
include("../../wavelab/winding_number.jl")

# Boussinesq driver


function boussinseq()

    # Parameters
    p = Parameter(0.4)

    # Profile
    s_I = 8
    s_R = s_I
    s_L = -s_I


    s = Infinity(s_I, s_R, s_L, A, Ak)

    # STABLAB structures

    # s, e, m, c = emcset(s, "front", [2,2], "default", A)
    # s, e, m, c = emcset(s, "front", [2,2], "reg_adj_polar", A)
    # s, e, m, c = emcset(s, "front", [2,2], "adj_reg_polar", A)
    s, e, m, c = emcset(s, "front", [2,2], "reg_reg_polar", A)
    # s, e, m, c = emcset(s, "front", [2,2], "reg_adj_compound", A, Ak)
    # s, e, m, c = emcset(s, "front", [2,2], "adj_reg_compound", A , Ak)


    # preimage

    points = 50
    preimage = 0.16 .+ 0.05 * exp.(2 * pi * 1im * LinRange(0, 0.5, points + (points - 1) * c.ksteps))


    # Compute Evans function

    halfw, _ = contour(c, s, p, m, e, preimage)
    w = [halfw[1 : end-1]; reverse(conj(halfw))]

    wnd = winding_number(w)
    println("Winding Number: ", wnd)
    plot(w/w[1])

end


boussinseq()