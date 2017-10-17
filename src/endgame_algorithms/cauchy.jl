export CauchyEndgame

struct CauchyEndgame <: AbstractEndgameAlgorithm
    samples_per_loop::Int
    loopclosed_tolerance::Float64
    # parameter for the first heuristic
    L::Float64
    # parameter for the second heuristic
    #β::Float64
    K::Float64
end

function CauchyEndgame(;
    samples_per_loop=8,
    loopclosed_tolerance=1e-6,
    L=0.75,
    #β=1e-8,
    K=0.5)
    CauchyEndgame(
        samples_per_loop,
        loopclosed_tolerance,
        L,
        #β,
        K)
end
