using Combinatorics

function wedgie(M)
    # Take wedgie product of columns of m

    n, k = size(M)

    combos = stack(collect(combinations(1:n, k)), dims = 1)
    
    # TODO:: I am pretty sure that sorting the rows like this is unnecessary. You would think combinations would generate in this order. Maybe check this out?
    combos = sortslices(combos, dims = 1)

    out = zeros(ComplexF64, size(combos, 1), 1)

    for j=axes(combos,1)
        out[j] = det(M[combos[j,:], :])

    end

    return out

end