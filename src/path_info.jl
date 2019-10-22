export path_info

"""
    path_info(tracker::CoreTracker, x, t₁ = 1.0, t₀ = 0.0)
    path_info(tracker::PathTracker, x)

Track a path using the given `tracker` and start value `x`. This returns a struct containing
detailed informations about the tracked path.

# Example
```julia
julia> @polyvar x y;
julia> f = [x^2 + y^2 + 1, x + y - 3];
julia> tracker, starts = coretracker_startsolutions(f);
julia> path_info(tracker, first(starts))
CTPathInfo:
 • # return code → success
 • # steps (✓/✗) → 11 ( 10 / 1 )
 • # factorizations → 31
 • # ldivs → 66
┌───┬───────┬────────┬───────┬─────────┬─────────┬──────┬─────────┬───────────┬──────┬─────────┐
│   │     s │     Δs │     ω │   |Δx₀| │     acc │    κ │       ψ │ limit_acc │  |x| │     |r| │
├───┼───────┼────────┼───────┼─────────┼─────────┼──────┼─────────┼───────────┼──────┼─────────┤
│ ✓ │     0 │  0.117 │     1 │       0 │ 3.4e-17 │  1.5 │ 2.2e-17 │   3.4e-17 │    1 │ 6.7e-17 │
│ ✓ │ 0.117 │  0.117 │ 0.866 │ 0.00063 │ 2.7e-17 │ 2.31 │ 1.7e-17 │   2.7e-17 │ 1.13 │ 5.6e-17 │
│ ✓ │ 0.235 │    0.2 │ 0.817 │ 0.00086 │ 3.7e-14 │ 2.31 │ 1.7e-17 │   2.7e-17 │ 1.29 │ 5.6e-17 │
│ ✗ │ 0.435 │  0.166 │ 0.437 │   0.021 │ 3.4e-09 │ 2.31 │ 1.7e-17 │   2.7e-17 │ 1.53 │ 5.6e-17 │
│ ✓ │ 0.435 │  0.105 │  1.08 │   0.041 │ 3.4e-09 │ 2.31 │ 1.7e-17 │   2.7e-17 │ 1.53 │ 5.6e-17 │
│ ✓ │ 0.539 │  0.105 │  1.08 │  0.0053 │   1e-15 │ 7.03 │ 6.7e-16 │     1e-15 │ 1.47 │ 2.2e-15 │
│ ✓ │ 0.644 │ 0.0704 │ 0.652 │   0.039 │ 3.4e-08 │ 7.03 │ 6.7e-16 │     1e-15 │ 1.24 │ 2.2e-15 │
│ ✓ │ 0.714 │ 0.0816 │ 0.408 │  0.0057 │ 1.7e-12 │ 7.03 │ 6.7e-16 │     1e-15 │ 1.33 │ 2.2e-15 │
│ ✓ │ 0.796 │ 0.0999 │ 0.357 │  0.0049 │ 3.3e-12 │ 7.03 │ 6.7e-16 │     1e-15 │ 1.59 │ 2.2e-15 │
│ ✓ │ 0.896 │  0.104 │ 0.683 │  0.0032 │ 4.1e-12 │ 7.03 │ 6.7e-16 │     1e-15 │ 1.94 │ 2.2e-15 │
│ ✓ │     1 │      0 │ 0.748 │ 0.00082 │ 1.5e-16 │ 3.86 │ 6.8e-17 │   1.5e-16 │ 2.24 │ 1.3e-15 │
└───┴───────┴────────┴───────┴─────────┴─────────┴──────┴─────────┴───────────┴──────┴─────────┘

julia> tracker = pathtracker(f);

julia> path_info(tracker, first(starts))
PTPathInfo:
 • # return code → success
 • # steps (✓/✗) → 19 ( 18 / 1 )
 • # factorizations → 56
 • # ldivs → 114
┌───┬───────┬───────┬─────────┬──────────┬──────────┬──────┬───────────┬─────────┬───────┬─────────┬──────┐
│   │     s │    Δs │    |ν̇| │    min_ν │    max_ν │    κ │ limit_acc │     acc │     ω │   |Δx₀| │  |x| │
├───┼───────┼───────┼─────────┼──────────┼──────────┼──────┼───────────┼─────────┼───────┼─────────┼──────┤
│ ✓ │     0 │ 0.126 │     NaN │      NaN │      NaN │  1.5 │   9.6e-18 │ 9.6e-18 │     1 │       0 │    1 │
│ ✓ │ 0.126 │ 0.126 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 5.4e-17 │ 0.874 │ 0.00098 │ 1.13 │
│ ✓ │ 0.251 │ 0.209 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 5.6e-14 │ 0.829 │ 0.00094 │ 1.29 │
│ ✓ │  0.46 │ 0.234 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │   3e-11 │ 0.733 │  0.0049 │ 1.51 │
│ ✓ │ 0.694 │  0.23 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 2.6e-10 │ 0.767 │   0.008 │ 1.54 │
│ ✓ │ 0.924 │ 0.158 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │   6e-08 │  1.16 │   0.028 │ 1.31 │
│ ✓ │  1.08 │ 0.154 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 6.3e-11 │ 0.598 │  0.0093 │ 1.21 │
│ ✓ │  1.24 │  0.25 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 8.9e-15 │ 0.395 │  0.0015 │  1.3 │
│ ✓ │  1.48 │  0.38 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │ 8.7e-14 │ 0.293 │  0.0022 │ 1.51 │
│ ✗ │  1.87 │ 0.135 │     NaN │      NaN │      NaN │  1.5 │   5.4e-17 │   8e-13 │ 0.686 │  0.0022 │ 1.78 │
│ ✓ │     2 │ 0.314 │    0.39 │    -0.26 │    -0.26 │  1.5 │   5.8e-17 │ 5.8e-17 │ 0.566 │ 8.8e-06 │ 1.85 │
│ ✓ │  2.31 │ 0.314 │    0.23 │    -0.17 │    -0.16 │ 3.64 │   1.2e-16 │ 1.2e-16 │ 0.491 │  0.0001 │ 1.97 │
│ ✓ │  2.63 │  1.24 │    0.15 │    -0.11 │     -0.1 │ 3.64 │   1.2e-16 │ 3.7e-10 │  0.84 │   3e-05 │ 2.05 │
│ ✓ │  3.86 │  1.71 │   0.028 │   -0.025 │   -0.023 │ 3.64 │   1.2e-16 │ 1.4e-12 │ 0.751 │   0.002 │ 2.19 │
│ ✓ │  5.58 │  2.84 │  0.0043 │  -0.0042 │  -0.0037 │ 3.71 │   7.9e-17 │ 7.9e-17 │ 0.275 │  0.0017 │ 2.23 │
│ ✓ │  8.41 │  5.13 │ 0.00024 │ -0.00024 │ -0.00021 │ 3.71 │   7.9e-17 │ 2.3e-14 │ 0.244 │  0.0012 │ 2.24 │
│ ✓ │  13.5 │  13.3 │ 1.4e-06 │ -1.4e-06 │ -1.3e-06 │ 3.71 │   7.9e-17 │ 1.1e-08 │ 0.245 │  0.0003 │ 2.24 │
│ ✓ │  26.8 │    55 │ 4.3e-10 │ -2.5e-12 │ -2.2e-12 │ 3.67 │   1.7e-16 │ 1.7e-16 │ 0.275 │ 9.8e-06 │ 2.24 │
│ ✓ │  81.8 │     0 │ 4.3e-10 │ -2.5e-12 │ -2.2e-12 │ 3.67 │   1.2e-16 │ 1.2e-16 │ 0.275 │ 1.8e-09 │ 2.24 │
└───┴───────┴───────┴─────────┴──────────┴──────────┴──────┴───────────┴─────────┴───────┴─────────┴──────┘
```
"""
function path_info end


struct CTPathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    Δx₀::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    cond::Vector{Float64}
    accuracy::Vector{Float64}
    residual::Vector{Float64}
    limit_accuracy::Vector{Float64}
    eval_err::Vector{Float64}
    # total info
    return_code::CoreTrackerStatus.states
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(tracker::CoreTracker, x₀, t₁ = 1.0, t₀ = 0.0)
    state = tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    Δx₀ = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    limit_accuracy = Float64[]
    residual = Float64[]
    eval_err = Float64[]

    init!(tracker, x₀, t₁, t₀)
    push!(s, state.s)
    push!(Δs, state.Δs)
    push!(ω, state.ω)
    push!(Δx₀, state.norm_Δx₀)
    push!(norm_x, maximum(abs, state.x))
    push!(cond, state.jacobian.cond[])
    push!(accuracy, state.accuracy)
    push!(limit_accuracy, state.limit_accuracy)
    push!(residual, maximum(tracker.corrector.abs_r))
    push!(eval_err, state.eval_err)

    first = true
    for _ in tracker
        push!(accepted_rejected, !state.last_step_failed)
        push!(s, state.s)
        push!(Δs, state.Δs)
        push!(ω, state.ω)
        push!(Δx₀, state.norm_Δx₀)
        push!(norm_x, maximum(abs, state.x))
        push!(cond, state.jacobian.cond[])
        push!(accuracy, state.accuracy)
        push!(limit_accuracy, state.limit_accuracy)
        push!(residual, maximum(tracker.corrector.abs_r))
        push!(eval_err, state.eval_err)
        first = false
    end
    push!(accepted_rejected, !state.last_step_failed)

    CTPathInfo(
        s,
        Δs,
        ω,
        Δx₀,
        accepted_rejected,
        norm_x,
        cond,
        accuracy,
        residual,
        limit_accuracy,
        eval_err,
        status(tracker),
        state.jacobian.factorizations[],
        state.jacobian.ldivs[],
    )
end

path_table(info::CTPathInfo) = path_table(stdout, info)
function path_table(io::IO, info::CTPathInfo)
    header = ["", "s", "Δs", "ω", "|Δx₀|", "acc", "κ", "ψ", "limit_acc", "|x|", "|r|"]
    h1 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 1 && data[i, 1] == :✗,
        crayon = PrettyTables.crayon"red",
    )
    h2 = PrettyTables.Highlighter(
        f = (data, i, j) -> j == 1 && data[i, 1] == :✓,
        crayon = PrettyTables.crayon"green",
    )
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    data = hcat(
        ✓✗,
        _sigdigits3.(info.s),
        _sigdigits3.(info.Δs),
        _sigdigits3.(info.ω),
        _sigdigits2.(info.Δx₀),
        _sigdigits2.(info.accuracy),
        _sigdigits3.(info.cond),
        _sigdigits2.(info.eval_err),
        _sigdigits2.(info.limit_accuracy),
        _sigdigits3.(info.norm_x),
        _sigdigits2.(info.residual),
    )
    PrettyTables.pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::CTPathInfo)
    println(io, "CTPathInfo:")
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


struct PTPathInfo
    # per step info
    s::Vector{Float64}
    Δs::Vector{Float64}
    ω::Vector{Float64}
    Δx₀::Vector{Float64}
    accepted_rejected::Vector{Bool}
    norm_x::Vector{Float64}
    cond::Vector{Float64}
    accuracy::Vector{Float64}
    limit_accuracy::Vector{Float64}
    norm_ν̇::Vector{Float64}
    min_ν::Vector{Float64}
    max_ν::Vector{Float64}
    # total info
    return_code::PathTrackerStatus.states
    n_factorizations::Int
    n_ldivs::Int
end

function path_info(tracker::PathTracker, x₀, t₀ = nothing)
    ct_state = tracker.core_tracker.state

    s = Float64[]
    Δs = Float64[]
    ω = Float64[]
    Δx₀ = Float64[]
    accepted_rejected = Bool[]
    norm_x = Float64[]
    cond = Float64[]
    accuracy = Float64[]
    limit_accuracy = Float64[]
    norm_ν̇ = Float64[]
    min_ν = Float64[]
    max_ν = Float64[]

    if t₀ !== nothing
        init!(tracker, x₀, 0.0, t₀)
    else
        init!(tracker, x₀)
    end

    push!(s, real(current_t(ct_state)))
    push!(Δs, real(current_Δt(ct_state)))
    push!(ω, ct_state.ω)
    push!(Δx₀, ct_state.norm_Δx₀)
    push!(norm_x, maximum(abs, ct_state.x))
    push!(cond, ct_state.jacobian.cond[])
    push!(accuracy, ct_state.accuracy)
    push!(limit_accuracy, ct_state.limit_accuracy)
    push!(norm_ν̇, maximum(abs, tracker.state.valuation.ν̇))
    push!(max_ν, maximum(tracker.state.valuation.ν))
    push!(min_ν, minimum(tracker.state.valuation.ν))

    first = true
    for _ in tracker
        push!(accepted_rejected, !ct_state.last_step_failed)
        push!(s, real(current_t(ct_state)))
        push!(Δs, real(current_Δt(ct_state)))
        push!(ω, ct_state.ω)
        push!(Δx₀, ct_state.norm_Δx₀)
        push!(norm_x, maximum(abs, ct_state.x))
        push!(cond, ct_state.jacobian.cond[])
        push!(accuracy, ct_state.accuracy)
        push!(limit_accuracy, ct_state.limit_accuracy)
        push!(norm_ν̇, maximum(abs, tracker.state.valuation.ν̇))
        push!(max_ν, maximum(tracker.state.valuation.ν))
        push!(min_ν, minimum(tracker.state.valuation.ν))
        first = false
    end
    push!(accepted_rejected, !ct_state.last_step_failed)

    PTPathInfo(
        s,
        Δs,
        ω,
        Δx₀,
        accepted_rejected,
        norm_x,
        cond,
        accuracy,
        limit_accuracy,
        norm_ν̇,
        min_ν,
        max_ν,
        status(tracker),
        ct_state.jacobian.factorizations[],
        ct_state.jacobian.ldivs[],
    )
end

path_table(info::PTPathInfo) = path_table(stdout, info)
function path_table(io::IO, info::PTPathInfo)
    header = [
        "",
        "s",
        "Δs",
        "|ν̇|",
        "min_ν",
        "max_ν",
        "κ",
        "limit_acc",
        "acc",
        "ω",
        "|Δx₀|",
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
    ✓✗ = map(v -> v ? :✓ : :✗, info.accepted_rejected)
    data = hcat(
        ✓✗,
        _sigdigits3.(info.s),
        _sigdigits3.(info.Δs),
        _sigdigits2.(info.norm_ν̇),
        _sigdigits2.(info.min_ν),
        _sigdigits2.(info.max_ν),
        _sigdigits3.(info.cond),
        _sigdigits2.(info.limit_accuracy),
        _sigdigits2.(info.accuracy),
        _sigdigits3.(info.ω),
        _sigdigits2.(info.Δx₀),
        _sigdigits3.(info.norm_x),
    )
    PrettyTables.pretty_table(io, data, header, crop = :none, highlighters = (h1, h2))
end

function Base.show(io::IO, info::PTPathInfo)
    println(io, "PTPathInfo:")
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

_sigdigits3(x) = Printf.@sprintf("%.3g", x)
_sigdigits2(x) = Printf.@sprintf("%.2g", x)
