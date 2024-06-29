function contour(c, s, p, m, e, pre_preimage)
    # Returns the Evans function output for the given input. The structures
    # c, s, p, m, and e are described in the STABLAB documentation. The 
    # input pre_preimage contains the contour points from which the Evans 
    # function will be evaluated. The sctructure c should contain a field 
    # lambda_steps and ksteps. The positive integer, c.ksteps indicates how
    # many Kato steps will be taken between points on which the Evans function
    # is initially evaluated. The positive integer c.lambda_steps indicates how
    # many additional points are specified between Kato steps in the contour,
    # pre_preimage, for optional evaluation of the Evans function if needed to 
    # obtain desired relative error, if specified. 
    #
    # Example: Suppose we want to evaluate the Evans function on a contour with
    # 3 points. In order to get accurate results we determine that we need to
    # take 2 Kato steps between each point. Then our preimage contour will have
    # entries [ 1 2 3 4 5 6 7] with entries 1, 4, and 7 corresponding to the 3
    # points we want to Evaluate the Evans function on and entries 2 and 3 intermediate
    # points on our contour between 1 and 4, and points 5 and 6 intermediate points on
    # the preimage contour between 4 and 7. If we want our Evans
    # function output to vary in relative error bewteen consecutive points by
    # less than some tolerance, c.tol, then we set c.refine = 'on'. Now the Evans
    # fucntion will be evaluated also on entries 2,3,5, and 6 if needed. Perhaps we
    # can achieve relative error only evaluating the Evans function on the
    # additional point 3. Then the Evans function is not evaluated on points
    # 2, 5, and 6. Suppose now that even after evaluating the Evans function on
    # all the preimage points we don't meet relative tolerance requirements.
    # Perhaps the region where we don't meet tolerance is only in one small
    # region. We don't want to slow the computation way down by computing the
    # Kato basis numerically on additional points on all of the contour. If we
    # specify c.lambda_steps = 1, for example, then between each Kato step, we
    # can compute the Evans function on an extra point if needed without
    # originally computing the Evans function on that point. So now our
    # preimage contour has 13 entries with the original points we compute on
    # residing in entries 1 7, and 13. The entries on which the Kato basis are
    # computed are 1,3,5,7,9,11,13. The entries on which the Kato basis can be
    # computed if needed and the Evans fucntion evaluated are 2,4,6,8,10,12.
    #
    # If desired, one can set c.check = 'on'. Then the Evans fucntion is first
    # evaluated on the last entry of the preimage contour to determine if the
    # output and the conjugate of the output are within specified relative
    # error, c.tol. This can save time computing the Evans fucntion on half
    # a contour approaching the origin if the Evans function is going to fail
    # at the origin. 
    #
    # If c.stats = 'on', then a waitbar showing computation 
    # progress is displayed. If c.stats = 'print', then computation progress is
    # printed to the command window instead of to a waitbar (use this option
    # for parallel computing)


    if cmp(c.stats, "on") == 0
        cstats = 1
    elseif cmp(c.stats, "print") == 0
        println("Finding the Kato basis")
    else
        cstats = 0
    end

    # Find the subset on which Kato basis is initially evaluated
    pre_index = 1:(c.lambda_steps + 1):length(pre_preimage)
    preimage = pre_preimage[pre_index]

    # Find the subset on which Evans function is initially evaluated
    lbasis, lproj = c.basisL(c.Lproj, c.L, preimage, s, p, c.LA, 1, c.epsl)
    rbasis, rproj = c.basisR(c.Rproj, c.R, preimage, s, p, c.RA, -1, c.epsr)
    #println("Lbasis: ", lbasis)
    #println("rbasis: ", rbasis)

    index = 1:(c.ksteps+1):length(preimage)
    lbasis2 = lbasis[:,:,index]
    rbasis2 = rbasis[:,:,index]
    preimage2 = preimage[index]

    out = zeros(ComplexF64, length(preimage2))

    # Makes sure the contour can be successfully computed close enough to the origin
    # to satisfy tolerance before computing everything

    #TODO:: Verify that this checking code actually works
    if c.check == "on"
        try
            #TODO:: There may be an error in the original MATLAB code on this line. Has preimage2[:,1] which doesn't work in Julia
            out[1] = c.evans(lbasis2[:,:,end], rbasis2[:,:,end], preimage2[1], s, p, m, e)
            near_origin = c.evans(lbasis2[:,:,end], rbasis2[:,:,end], preimage2[:,end], s, p, m, e)
            out[end] = near_origin / out[1]

            if abs(conj(out[end]) - out[end]) / abs(out[end]) > c.tol
                println(out[1])
                println(out[end])
                println(abs(conj(out[end]) - out[end]) / abs(out[end]))
                error("The Evans function does not satisfy tolerance at the endpoint of the contour")
            end
        catch e
            error("The Evans function failed to compute at the end point")
        end
    end


    # Compute the Evans function on the initial contour
    if cstats == 1
        # TODO:: MATLAB code uses a waitbar here. Decide if we want to attempt something similar here
        # TODO:: MATLAB code begins timing things here. Decide how to best implement a timing mechanism in Julia

        for j = 1:length(index)
            out[j] = c.evans(lbasis2[:, :, j], rbasis[:, :, j], preimage2[:, j], s, p, m, e)
        end

    else
        if c.stats == "print"
            println("Computing the Evans function on the first set of points")
        end
        if c.debug == "on"
            for j = 1:length(index)
                out[j] = c.evans(lbasis2[:, :, j], rbasis2[:, :, j], preimage2[j], s, p, m, e)
            end
        else
            # TODO:: MATLAB code uses a parallelized loop here. Looks like there are ways to do this Julia but probably want to research more about how it works and how stable the package is across platforms https://discourse.julialang.org/t/trying-to-write-a-parallel-for-loop-in-julia/40862/3
            for j = 1:length(index)
                out[j] = c.evans(lbasis2[:, :, j], rbasis2[:, :, j], preimage2[j], s, p, m, e)
            end
        end

    end

    # Check if relative error tolerance has been specified
    if c.refine != "on"
        return out, preimage2
    end

    # Refine the mesh on which the Evans function is computed until requested tolerance is achieved
    # using the Kato steps as needed

    # TODO:: Finish the rest of the function



    # TODO:: Modify return statement when function is complete
    return 0, 0
    #return out, preimage2

end



