include("A.jl")

include("../../wavelab/structs.jl")
include("../../wavelab/emcset.jl")
include("../../wavelab/semicirc2.jl")
include("../../wavelab/contour.jl")

# TODO::Is there a Beep off equivalent in Julia? Do we care enough to do it?


# Parameters
p = Parameter(1, 0, "off")

#Numerical Infinity
s = Infinity(12, 12, -12, nothing, nothing)

# set wavelab structures to local default values
s, e, m, c = emcset(s, "front", [1, 1], "default", A)   # default for Burgers is reg_reg_polar

# preimage contour
circpnts = 20
imagpnts = 20
innerpnts = 5
r = 10
spread = 4
zerodist = 1e-2

preimage = semicirc2(circpnts, imagpnts, innerpnts, c.ksteps, r, spread, zerodist, c.lambda_steps)



#Compute the Evans function
halfw = contour(c, s, p, m, e, preimage)



