# Boussinesq driver


function boussinseq()

    # Parameters
    p.S = 0.4

    # Profile
    s.I = 8
    s.R = s.I
    s.L = -s.I

    # STABLAB structures

    # s, e, m, c = emcset(s, "front", [2,2], 'default')
    # s, e, m, c = emcset(s, "front", [2,2], 'reg_adj_polar')
    # s, e, m, c = emcset(s, "front", [2,2], 'adj_reg_polar')
    # s, e, m, c = emcset(s, "front", [2,2], 'reg_reg_polar')
    s, e, m, c = emcset(s, "front", [2,2], 'reg_adj_compound')
    # s, e, m, c = emcset(s, "front", [2,2], 'adj_reg_compound')


    # preimage

    points = 50
    preimage = 0.16 + 0.05 * exp(2 * pi * 1im * LinRange(0, 0.5, points + (points - 1) * c.ksteps))


    # Compute Evans function

    halfw = contour(c.s, p, m, e, preimage)
    w = [halfw[1 : end-1]; fliplr(conj(halfw))]

    wnd = winding_number(w)
    plot(w/w[1])
    println("Winding Number: " + wnd)

end


boussinseq()