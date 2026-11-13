# Stdlib-only loader for benchmarks/systems/*.toml, shared by the harness scripts.
# `build_system` takes the `Variable` and `System` constructors to use (HC.jl's ModelKit).
using TOML

const SYSROOT = normpath(joinpath(@__DIR__, "..", "systems"))

function _terms!(acc, ex, sign = 1)
    if ex isa Expr && ex.head == :call && ex.args[1] === :+
        foreach(a -> _terms!(acc, a, sign), ex.args[2:end])
    elseif ex isa Expr && ex.head == :call && ex.args[1] === :- && length(ex.args) == 3
        _terms!(acc, ex.args[2], sign); _terms!(acc, ex.args[3], -sign)
    else
        push!(acc, sign == 1 ? ex : :(-$ex))
    end
    acc
end

"""
    load_system(path) -> (; name, vars, params, eqs)

`path` is relative to benchmarks/systems. For `kind = "coefficients"` systems without
parameters every top-level term gets a fresh coefficient parameter.
"""
function load_system(path)
    d = TOML.parsefile(joinpath(SYSROOT, path))
    vars = Symbol.(d["variables"]); params = Symbol.(d["parameters"])
    eqs = Any[Meta.parse(s) for s in d["equations"]]
    if isempty(params) && d["parametrization"]["kind"] == "coefficients"
        # coefficient names c1, c2, … (c_1, … if a variable is already called like that)
        c = any(v -> occursin(r"^c[0-9]+$", String(v)), vars) ? "c_" : "c"
        k = 0
        eqs = map(eqs) do e
            Expr(:call, :+, [(k += 1; :($(Symbol(c, k)) * $t)) for t in _terms!(Any[], e)]...)
        end
        params = [Symbol(c, i) for i = 1:k]
    end
    (; name = replace(basename(path), ".toml" => ""), vars, params, eqs)
end

"""
    instance_parameters(path) -> Vector{ComplexF64}

The parameters of the system's stored instance: `target_parameters`, or all coefficients 1
for a coefficient family (the stored equations themselves).
"""
function instance_parameters(path)
    d = TOML.parsefile(joinpath(SYSROOT, path))
    haskey(d, "target_parameters") && return parse.(ComplexF64, d["target_parameters"])
    ones(ComplexF64, length(load_system(path).params))
end

_ev(x::Number, env) = x
_ev(x::Symbol, env) = x === :im ? im : env[x]
function _ev(ex::Expr, env)
    op = ex.args[1]
    a = [_ev(y, env) for y in ex.args[2:end]]
    op === :^ && return a[1]^Int(a[2])
    op === :+ && return reduce(+, a)
    op === :* && return reduce(*, a)
    op === :- && return length(a) == 1 ? -a[1] : a[1] - a[2]
    op === :/ && return a[1] / a[2]
    error("unsupported op $op")
end

function build_system(sys, Variable, System)
    xs = [Variable(v) for v in sys.vars]; ps = [Variable(v) for v in sys.params]
    env = Dict{Symbol,Any}(zip(sys.vars, xs)); merge!(env, Dict(zip(sys.params, ps)))
    System([_ev(e, env) for e in sys.eqs]; variables = xs, parameters = ps)
end

# Systems in the harness: path relative to benchmarks/systems.
const HARNESS_SYSTEMS = [
    "tier_b/hyperbola_path_jumping.toml",
    "tier_b/coinciding_paths_no_branch_point.toml",
    "tier_a/five_point_relative_pose.toml",
    "tier_a/six_revolute_ik.toml",
    "tier_a/stewart_gough_fk.toml",
    "tier_b/eco8.toml",
    "tier_b/noon5.toml",
    "tier_b/cyclic7.toml",
    "tier_b/katsura10.toml",
    "hcjl/cyclooctane.toml",
    "tier_a/steiner_conics.toml",
    # from HC.jl issues (random → random here; → instance in benchmarks/harness/reference.jl)
    "tier_b/zero_solution_constant_path_issue647.toml",
    "tier_a/biquadratic_factorization_issue649.toml",
    "tier_c/bicritical_cycle_issue502.toml",
    "tier_c/rc_circuit_badly_scaled_issue414.toml",
    "tier_c/jacobian_3torsion_genus3_issue574.toml",
    # hard: full solution sets from benchmarks/harness/solution_sets.jl
    "hcjl/tritangents.toml",
    "hcjl/fano_quintic_lines.toml",
    "tier_a/alt_nine_point_four_bar.toml",
]
