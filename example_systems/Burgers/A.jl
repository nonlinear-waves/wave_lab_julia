function A(x, lambda, s, p)
    # Evans matrix for Burgers system in unintegrated coordinates

    a = 0.5 * (p.ul - p.ur)
    cc = 0.5 * (p.ul + p.ur)  # Wave speed 
    u = cc - a * tanh(a * x / 2) # Profile
    uder = (-a^2 / 2) * sech(a * x / 2)^2 # profile derivative

    if cmp(p.integrated, "on") == 0
        out = [0 1 
               lambda u-cc]
    else
        out = [0 1
               lambda+uder u-cc]
    end

    
    return out
end