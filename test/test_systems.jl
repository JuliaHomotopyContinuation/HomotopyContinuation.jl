function cyclic(n)
    @var z[1:n]
    eqs = [sum(prod(z[(k-1)%n+1] for k = j:j+m) for j = 1:n) for m = 0:(n-2)]
    push!(eqs, prod(z) - 1)
    System(eqs, z)
end

function bacillus()
    @var w w2 w2v v w2v2 vP sigmaB w2sigmaB vPp phos
    poly = [
        (-1 * 0.7 * w + -2 * 3600.0 * (w^2 / 2) + 2 * 18.0 * w2) * (0.2 + sigmaB) +
        4.0 * 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2 +
        3600.0 * (w^2 / 2) +
        -1 * 18.0 * w2 +
        -1 * 3600.0 * w2 * v +
        18.0w2v +
        36.0w2v +
        -1 * 3600.0 * w2 * sigmaB +
        18.0w2sigmaB,
        -1 * 0.7 * w2v +
        3600.0 * w2 * v +
        -1 * 18.0 * w2v +
        -1 * 3600.0 * w2v * v +
        18.0w2v2 +
        -1 * 36.0 * w2v +
        36.0w2v2 +
        1800.0 * w2sigmaB * v +
        -1 * 1800.0 * w2v * sigmaB,
        (
            -1 * 0.7 * v +
            -1 * 3600.0 * w2 * v +
            18.0w2v +
            -1 * 3600.0 * w2v * v +
            18.0w2v2 +
            -1 * 1800.0 * w2sigmaB * v +
            1800.0 * w2v * sigmaB +
            180.0vPp
        ) * (0.2 + sigmaB) + 4.5 * 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2v2 + 3600.0 * w2v * v + -1 * 18.0 * w2v2 + -1 * 36.0 * w2v2,
        -1 * 0.7 * vP + 36.0w2v + 36.0w2v2 + -1 * 3600.0 * vP * phos + 18.0vPp,
        (
            -1 * 0.7 * sigmaB +
            -1 * 3600.0 * w2 * sigmaB +
            18.0w2sigmaB +
            1800.0 * w2sigmaB * v +
            -1 * 1800.0 * w2v * sigmaB
        ) * (0.2 + sigmaB) + 0.4 * (1 + 30.0sigmaB),
        -1 * 0.7 * w2sigmaB +
        3600.0 * w2 * sigmaB +
        -1 * 18.0 * w2sigmaB +
        -1 * 1800.0 * w2sigmaB * v +
        1800.0 * w2v * sigmaB,
        -1 * 0.7 * vPp + 3600.0 * vP * phos + -1 * 18.0 * vPp + -1 * 180.0 * vPp,
        (phos + vPp) - 2.0,
    ]
    System(poly, [w, w2, w2v, v, w2v2, vP, sigmaB, w2sigmaB, vPp, phos])
end
