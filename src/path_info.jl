export path_info

struct PathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    Δx₀::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    accuracy::Vector{Float64}
    μ::Vector{Float64}
    τ::Vector{Float64}
    Δx_t::Vector{Float64}
    Δx̂x::Vector{Float64}
    high_prec::Vector{Bool}
    # total info
    return_code::TrackerCondition.conditions
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(
    tracker::Tracker,
    x₀,
    t₁ = 0.0,
    t₀ = 1.0;
    debug::Bool = false,
)
    state = tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    Δx₀ = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    μ = Float64[]
    τ = Float64[]
    Δx_t = Float64[]
    Δx̂x = Float64[]
    high_prec = Bool[]

    p = order(tracker.predictor)

    init!(tracker, x₀, t₁, t₀)
    while is_tracking(tracker.state.condition)
        push!(s, real(state.t))
        push!(Δs, real(state.Δt))
        e = tracker.state.norm(local_error(tracker.predictor))
        push!(Δx_t, e * abs(state.Δt)^p)
        push!(τ, tracker.options.β_τ * trust_region(tracker.predictor))
        push!(ω, state.ω)
        push!(accuracy, state.accuracy)
        push!(μ, state.μ)
        push!(high_prec, state.high_prec_residual)
        step!(tracker, debug)

        push!(accepted_rejected, !state.last_step_failed)
        push!(Δx₀, state.norm_Δx₀)
        push!(norm_x, maximum(abs, state.x))
        push!(Δx̂x, state.norm(state.x .- state.x̂))
    end

    PathInfo(
        s,
        Δs,
        ω,
        Δx₀,
        accepted_rejected,
        norm_x,
        accuracy,
        μ,
        τ,
        Δx_t,
        Δx̂x,
        high_prec,
        tracker.state.condition,
        state.jacobian.factorizations[],
        state.jacobian.ldivs[],
    )
end

path_table(info::PathInfo) = path_table(stdout, info)
function path_table(io::IO, info::PathInfo)
    header = [
        "",
        "s",
        "Δs",
        "ω",
        "|Δx₀|",
        "h₀",
        "acc",
        "μ",
        "τ",
        "Δx_t",
        "Δpred ",
        "|x|",
    ]
    h1 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 1 && data[i, 1] == :✗,
        crayon = PrettyTables.crayon"red",
    )
    h2 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 1 && data[i, 1] == :✓,
        crayon = PrettyTables.crayon"green",
    )
    h3 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 3 && data[i, 3] ≈ -data[i, 9],
        crayon = PrettyTables.crayon"blue",
    )
    h4 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 8 && info.high_prec[i],
        crayon = PrettyTables.crayon"blue",
    )

    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    data = hcat(
        ✓✗,
        info.s,
        info.Δs,
        info.ω,
        info.Δx₀,
        info.ω .* info.Δx₀,
        info.accuracy,
        info.μ,
        info.τ,
        info.Δx_t,
        info.Δx̂x,
        info.norm_x,
    )

    ft1 = PrettyTables.ft_printf("%5.2g", [4, 5, 6, 7, 8, 10, 11, 12])
    ft2 = PrettyTables.ft_printf("%5.2g", [3, 9])
    ft4 = PrettyTables.ft_printf("%3.3g", [2, 4])
    ft = merge(ft1, ft2, ft4)

    PrettyTables.pretty_table(
        io,
        data,
        header,
        crop = :none,
        formatter = ft,
        highlighters = (h1, h2, h3, h4),
    )
end

function Base.show(io::IO, info::PathInfo)
    println(io, "PathInfo:")
    println(io, " • # return code → ", info.return_code)
    println(
        io,
        " • # steps (✓/✗) → ",
        length(info.s),
        " ( ",
        count(info.accepted_rejected),
        " / ",
        count(!, info.accepted_rejected),
        " )",
    )
    println(io, " • # factorizations → ", info.n_factorizations)
    println(io, " • # ldivs → ", info.n_ldivs)
    path_table(io, info)
end
