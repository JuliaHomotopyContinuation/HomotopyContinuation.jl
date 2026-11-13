# Validate the benchmark system collection.
#
#     julia benchmarks/systems/validate.jl
#
# Needs only the Julia standard library (TOML). For every `<tier>/<name>.toml`
# it checks the schema described in README.md, that all equations parse with
# `Meta.parse`, use only + - * / ^ sqrt and declared names, that systems are
# square (unless `square = false`), and that every equation can be evaluated at
# a random complex point.

using TOML
using Random

const ROOT = @__DIR__
const TIERS = Dict("tier_a" => "A", "tier_b" => "B", "tier_c" => "C", "tier_d" => "D", "hcjl" => "HC")
const KINDS = ("explicit", "coefficients", "pointer", "fixed")
const OPS = (:+, :-, :*, :/, :^, :sqrt)

symbols!(acc, x::Symbol) = push!(acc, x)
symbols!(acc, ::Number) = acc
function symbols!(acc, ex::Expr)
    ex.head == :call || error("unsupported expression head `$(ex.head)` in `$ex`")
    ex.args[1] in OPS || error("unsupported operation `$(ex.args[1])` in `$ex`")
    foreach(a -> symbols!(acc, a), ex.args[2:end])
    acc
end

evaluate(x::Number, env) = x
evaluate(x::Symbol, env) = x === :im ? im : env[x]
function evaluate(ex::Expr, env)
    op = ex.args[1]
    args = [evaluate(a, env) for a in ex.args[2:end]]
    op === :^ && return args[1]^Int(args[2])
    op === :sqrt && return sqrt(complex(args[1]))
    return getfield(Base, op)(args...)
end

function parse_equation(s::AbstractString)
    ex = Meta.parse(s; raise = true)
    ex isa Expr && ex.head in (:incomplete, :error) && error("cannot parse `$s`")
    ex
end

function check_exponents(ex)
    ex isa Expr || return
    if ex.args[1] === :^
        e = ex.args[3]
        e isa Integer && e >= 0 || error("exponent must be a nonnegative integer literal in `$ex`")
    end
    foreach(check_exponents, ex.args[2:end])
end

function check_equations(eqs, vars, params, label, errs)
    names = Set(Symbol.([vars; params]))
    env = Dict{Symbol,ComplexF64}(n => randn(ComplexF64) for n in names)
    used = Set{Symbol}()
    nerrs = length(errs)
    for (i, s) in enumerate(eqs)
        try
            ex = parse_equation(s)
            check_exponents(ex)
            syms = symbols!(Set{Symbol}(), ex)
            unknown = setdiff(syms, names, (:im,))
            isempty(unknown) || error("undeclared names $(collect(unknown))")
            union!(used, syms)
            v = evaluate(ex, env)
            isfinite(v) || error("evaluates to $v at a random point")
        catch e
            push!(errs, "$label[$i]: $(sprint(showerror, e))")
        end
    end
    if length(errs) == nerrs # occurrence check is only meaningful if all equations parsed
        for v in vars
            Symbol(v) in used || push!(errs, "$label: variable `$v` does not occur in any equation")
        end
    end
    used
end

function check_values(vals, params, key, errs)
    length(vals) == length(params) ||
        push!(errs, "$key has $(length(vals)) entries, expected $(length(params))")
    for v in vals
        try
            ex = parse_equation(string(v))
            isempty(setdiff(symbols!(Set{Symbol}(), ex), (:im,))) || error("not a number")
            evaluate(ex, Dict{Symbol,ComplexF64}())
        catch e
            push!(errs, "$key: cannot evaluate `$v` as a number ($(sprint(showerror, e)))")
        end
    end
end

isname(s) = s isa String && Base.isidentifier(s) && s != "im"

function validate(path)
    errs = String[]
    d = try
        TOML.parsefile(path)
    catch e
        return ["TOML parse error: $(sprint(showerror, e))"], nothing
    end
    tierdir = basename(dirname(path))
    for key in ("name", "tier", "description", "tags", "variables", "parameters", "equations",
                "source", "solutions", "parametrization")
        haskey(d, key) || push!(errs, "missing key `$key`")
    end
    isempty(errs) || return errs, d

    d["name"] == splitext(basename(path))[1] || push!(errs, "name does not match file name")
    get(TIERS, tierdir, nothing) == d["tier"] || push!(errs, "tier `$(d["tier"])` does not match directory `$tierdir`")
    isempty(d["tags"]) && push!(errs, "tags must not be empty")
    for k in ("citation", "url")
        haskey(d["source"], k) || push!(errs, "missing source.$k")
    end

    vars, params, eqs = d["variables"], d["parameters"], d["equations"]
    for n in [vars; params]
        isname(n) || push!(errs, "invalid name `$n`")
    end
    allunique([vars; params]) || push!(errs, "variable/parameter names are not unique")

    kind = get(d["parametrization"], "kind", "")
    kind in KINDS || push!(errs, "parametrization.kind must be one of $KINDS")

    sol = d["solutions"]
    haskey(sol, "count_source") || push!(errs, "missing solutions.count_source")
    if get(sol, "zero_dimensional", true) == false
        get(d, "square", true) == false || push!(errs, "zero_dimensional = false requires square = false")
    elseif haskey(sol, "generic_nonsingular")
        sol["generic_nonsingular"] isa Integer && sol["generic_nonsingular"] > 0 ||
            push!(errs, "solutions.generic_nonsingular must be a positive integer")
    elseif !haskey(sol, "generic_nonsingular_formula") && get(sol, "count_known", true)
        push!(errs, "missing solutions.generic_nonsingular")
    end

    if kind == "pointer"
        haskey(d, "external") || push!(errs, "pointer systems need an [external] table")
    else
        length(eqs) == length(vars) || get(d, "square", true) == false ||
            push!(errs, "$(length(eqs)) equations for $(length(vars)) variables (not square)")
        check_equations(eqs, vars, params, "equations", errs)
    end

    if kind == "coefficients"
        isempty(params) || push!(errs, "coefficient-parametrized systems keep `parameters` empty")
        haskey(d["parametrization"], "description") || push!(errs, "missing parametrization.description")
        p = d["parametrization"]
        if haskey(p, "generic_equations")
            gp = get(p, "generic_parameters", String[])
            length(p["generic_equations"]) == length(vars) || push!(errs, "generic_equations not square")
            allunique([vars; gp]) || push!(errs, "generic_parameters clash with variables")
            used = check_equations(p["generic_equations"], vars, gp, "generic_equations", errs)
            for q in gp
                Symbol(q) in used || push!(errs, "generic parameter `$q` unused")
            end
        end
    end

    if haskey(d, "variable_groups")
        g = reduce(vcat, d["variable_groups"]; init = String[])
        sort(g) == sort(vars) || push!(errs, "variable_groups is not a partition of variables")
    end
    for key in ("example_parameters", "start_parameters", "target_parameters", "alternative_target_parameters")
        haskey(d, key) && check_values(d[key], params, key, errs)
    end
    for key in ("start_solution", "target_solution")
        haskey(d, key) && check_values(d[key], vars, key, errs)
    end
    if d["tier"] == "C"
        # a coefficient family's target is the stored instance itself (all coefficients 1)
        haskey(d, "target_parameters") || (isempty(params) && kind == "coefficients") ||
            push!(errs, "tier C systems need target_parameters")
        haskey(sol, "target_endpoints") || push!(errs, "tier C systems need solutions.target_endpoints")
    end
    return errs, d
end

function main()
    Random.seed!(0x1234)
    files = sort([joinpath(ROOT, t, f) for t in keys(TIERS) if isdir(joinpath(ROOT, t))
                  for f in readdir(joinpath(ROOT, t)) if endswith(f, ".toml")])
    nbad = 0
    for f in files
        errs, d = validate(f)
        rel = relpath(f, ROOT)
        if isempty(errs)
            n = get(d["solutions"], "generic_nonsingular", get(d["solutions"], "generic_nonsingular_formula", "?"))
            println("ok    $rel  ($(length(d["variables"])) vars, $(length(d["parameters"])) params, $n solutions)")
        else
            nbad += 1
            println("FAIL  $rel")
            foreach(e -> println("        ", e), errs)
        end
    end
    println("\n$(length(files) - nbad)/$(length(files)) systems valid")
    exit(nbad == 0 ? 0 : 1)
end

main()
