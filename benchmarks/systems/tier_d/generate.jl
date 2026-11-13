# Generate the tier D ("large") benchmark systems.
#
#     julia benchmarks/systems/tier_d/generate.jl
#
# Needs only the Julia standard library. Writes every `tier_d/<name>.toml`
# (overwriting existing files). The size of each family is set in the call at
# the bottom of the file; the generator functions take the size as argument,
# so other members of the families (katsura-n, cyclic-n, dense n/d, ...) can be
# produced the same way.

const OUT = @__DIR__

# ---------------------------------------------------------------- TOML output

q(s) = "\"" * replace(string(s), "\\" => "\\\\", "\"" => "\\\"") * "\""

function names_block(io, key, names)
    if isempty(names)
        println(io, "$key = []")
        return
    end
    println(io, "$key = [")
    for chunk in Iterators.partition(names, 8)
        println(io, "  ", join(q.(chunk), ", "), ",")
    end
    println(io, "]")
end

function lines_block(io, key, lines)
    println(io, "$key = [")
    foreach(l -> println(io, "  ", q(l), ","), lines)
    println(io, "]")
end

function table(io, name, pairs)
    println(io, "\n[$name]")
    for (k, v) in pairs
        if v isa AbstractVector
            names_block(io, k, v)
        elseif v isa AbstractString
            println(io, "$k = ", q(v))
        else
            println(io, "$k = ", v)
        end
    end
end

function write_system(s)
    path = joinpath(OUT, s.name * ".toml")
    open(path, "w") do io
        println(io, "name = ", q(s.name))
        println(io, "tier = \"D\"")
        println(io, "description = ", q(s.description))
        println(io, "tags = [")
        foreach(t -> println(io, "  ", q(t), ","), s.tags)
        println(io, "]")
        names_block(io, "variables", s.variables)
        names_block(io, "parameters", s.parameters)
        lines_block(io, "equations", s.equations)
        haskey(s, :example_parameters) && names_block(io, "example_parameters", s.example_parameters)
        haskey(s, :generator) && println(io, "generator = ", q(s.generator))
        table(io, "source", s.source)
        table(io, "solutions", s.solutions)
        table(io, "parametrization", s.parametrization)
        table(io, "verification", s.verification)
    end
    println("wrote $(relpath(path, dirname(OUT)))  ($(filesize(path)) bytes)")
end

const GENERATOR = "benchmarks/systems/tier_d/generate.jl"

# ------------------------------------------------------------------- helpers

# Monomials of degree <= d in n variables, ordered by total degree descending,
# then lexicographically descending in the exponent vector (x1 before x2 ...).
# The constant monomial is last.
function monomials(n, d)
    ms = Vector{Vector{Int}}()
    rec(prefix, left) = length(prefix) == n - 1 ? push!(ms, [prefix; left]) :
        for k in left:-1:0; rec([prefix; k], left - k); end
    for deg in d:-1:0
        n == 1 ? push!(ms, [deg]) : rec(Int[], deg)
    end
    ms
end

function monomial_string(vars, e)
    parts = String[]
    for (v, k) in zip(vars, e)
        k == 0 && continue
        push!(parts, k == 1 ? v : "$v^$k")
    end
    join(parts, "*")
end

function term_string(coef, mono)
    isempty(mono) && return coef
    isempty(coef) && return mono
    coef * "*" * mono
end

# Join signed terms (coefficient, rest) into an expression string.
function signed_sum(terms)
    io = IOBuffer()
    for (i, (c, m)) in enumerate(terms)
        neg = c < 0
        a = abs(c)
        body = m == "" ? string(a) : (a == 1 ? m : "$a*$m")
        if i == 1
            print(io, neg ? "-" : "", body)
        else
            print(io, neg ? " - " : " + ", body)
        end
    end
    String(take!(io))
end

const COEF_DESCRIPTION = "All coefficients are parameters: the generic family multiplies every monomial term of every equation (in the expanded form of the stored equations, including constants) by an independent complex parameter. The stored equations are one member of that family."

# ------------------------------------------------------------------ katsura

function katsura(n)
    vars = ["x$i" for i in 0:n]
    eqs = String[]
    for m in 0:n-1
        c = Dict{Tuple{Int,Int},Int}()
        for l in -n:n
            abs(m - l) <= n || continue
            i, j = minmax(abs(l), abs(m - l))
            c[(i, j)] = get(c, (i, j), 0) + 1
        end
        terms = [(c[k], k[1] == k[2] ? "x$(k[1])^2" : "x$(k[1])*x$(k[2])") for k in sort(collect(keys(c)))]
        push!(terms, (-1, "x$m"))
        push!(eqs, signed_sum(terms))
    end
    push!(eqs, signed_sum([(1, "x0"); [(2, "x$i") for i in 1:n]; (-1, "")]))
    vars, eqs
end

function katsura_system(n; verification)
    vars, eqs = katsura(n)
    (name = "katsura$n",
     description = "Katsura-$n (magnetism / mean-field spin glass), the same family as tier_b/katsura10 and tier_b/katsura12 at a larger size. n+1 = $(n+1) unknowns x0..x$n (with x_{-i} = x_i): for m = 0..$(n-1): sum_{l = -n..n, |m-l| <= n} x_l x_{m-l} - x_m = 0, plus x0 + 2 sum_{i=1..n} x_i - 1 = 0. $(n+1) quadrics with about n/2 terms each; Bezout number 2^$n = $(2^n).",
     tags = ["large", "stress", "scaling", "sparse", "coefficient-parameter", "generator"],
     variables = vars, parameters = String[], equations = eqs,
     generator = "katsura(n) in $GENERATOR: variables x0..xn; for m in 0..n-1: sum(x[abs(l)]*x[abs(m-l)] for l in -n..n if abs(m-l) <= n) - x[m]; then x[0] + 2*sum(x[1:n]) - 1. Identical to PolynomialTestSystems.jl katsura(n) and to tier_b/katsura12 (the generator reproduces that file's equations for n = 12).",
     source = ["citation" => "S. Katsura, 'Users posing problems to PoSSo', PoSSo Newsletter 2, 1994; S. Katsura et al., Cell Biophysics 11 (1987) 309-319. Standard form from PolynomialTestSystems.jl katsura(n). Standard scaling benchmark (HC.jl, Bertini 2 comparison, algpath-bench, Telen et al. go up to katsura-20; see tier_b/katsura12).",
               "url" => "https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/e04087ef08cf91ffafd88546c4d4ccb25613a3c7/src/systems.jl#L101-L109",
               "file" => "PolynomialTestSystems.jl@e04087ef08cf91ffafd88546c4d4ccb25613a3c7:src/systems.jl (katsura); generated by $GENERATOR"],
     solutions = ["generic_nonsingular" => 2^n,
                  "count_source" => "Bezout number 2^n = $(2^n); the Katsura systems attain it with all solutions finite and nonsingular (PHCpack Demo/katsura5..katsura10 list 2^n solutions; tier_b/katsura12 verified 4096 with HC.jl). As for katsura12, the torus mixed volume is smaller than 2^n because some solutions have zero coordinates, so a purely toric polyhedral start system misses solutions."],
     parametrization = ["kind" => "coefficients", "description" => COEF_DESCRIPTION * " Generic count of the family = 2^n."],
     verification = verification)
end

# ------------------------------------------------------------------- cyclic

function cyclic(n)
    vars = ["z$i" for i in 0:n-1]
    eqs = String[]
    for m in 0:n-2
        push!(eqs, join([join(["z$(mod(j + k, n))" for k in 0:m], "*") for j in 0:n-1], " + "))
    end
    push!(eqs, join(vars, "*") * " - 1")
    vars, eqs
end

function cyclic_system(n; mv, verification)
    vars, eqs = cyclic(n)
    (name = "cyclic$n",
     description = "Cyclic $n-roots: elementary cyclic sums sum_j prod_{k=0..m} z_{(j+k) mod n} = 0 for m = 0..$(n-2) and z0*...*z$(n-1) = 1. $n equations with $n terms each and degrees up to $n (total degree $(n)! = $(factorial(n))).",
     tags = ["large", "stress", "coefficient-parameter", "symmetric", "polyhedral", "high-degree"],
     variables = vars, parameters = String[], equations = eqs,
     generator = "cyclic(n) in $GENERATOR. Same system as tier_b/cyclic7 and HC.jl test/test_systems.jl cyclic(n).",
     source = ["citation" => "G. Bjorck, R. Froberg, 'A faster way to count the solutions of inhomogeneous systems of algebraic equations, with applications to cyclic n-roots', J. Symb. Comp. 12 (1991) 329-336. PHCpack Demo/cyclic10. Isolated solutions of the specific system: Y. Dai, S. Kim, M. Kojima, 'Computing all nonsingular solutions of cyclic-n polynomial using polyhedral homotopy continuation methods', J. Comput. Appl. Math. 152 (2003) 83-97.",
               "url" => "https://github.com/janverschelde/PHCpack/blob/ed77cc0234aec45099d1d587a90a1617d5bfdacd/src/Demo/cyclic$n",
               "file" => "PHCpack@ed77cc0234aec45099d1d587a90a1617d5bfdacd:src/Demo/cyclic$n (same equations, PHCpack lists them in the same order)"],
     solutions = ["generic_nonsingular" => mv,
                  "count_source" => "Coefficient family: Bernstein's theorem, generic count = mixed volume = $mv (PHCpack Demo/cyclic$n: 'BKK bound : 35940'). All solutions lie in the torus because z0*...*z$(n-1) = 1, so torus and affine counts agree. The specific cyclic-$n system (all coefficients 1) has 34940 isolated nonsingular solutions (Dai-Kim-Kojima 2003), so 1000 paths of a polyhedral homotopy diverge or end singular there."],
     parametrization = ["kind" => "coefficients", "description" => COEF_DESCRIPTION * " Generic count of the family = mixed volume = $mv."],
     verification = verification)
end

# ---------------------------------------------------------------------- eco

function eco(n)
    vars = ["x$i" for i in 1:n]
    eqs = String[]
    for k in 1:n-1
        inner = ["x$k"; ["x$i*x$(i+k)" for i in 1:n-k-1]]
        push!(eqs, "(" * join(inner, " + ") * ")*x$n - $k")
    end
    push!(eqs, join(vars[1:n-1], " + ") * " + 1")
    vars, eqs
end

function eco_system(n; count, verification)
    vars, eqs = eco(n)
    (name = "eco$n",
     description = "$n-dimensional economics problem (Morgan), the same family as tier_b/eco8 at a larger size: (x_k + sum_{i=1..n-k-1} x_i x_{i+k}) x_n = k for k = 1..n-1, and x_1 + ... + x_{n-1} + 1 = 0. Total degree 2*3^(n-2) = $(2*3^(n-2)) against $count solutions, so a total-degree homotopy has almost only diverging paths.",
     tags = ["large", "stress", "coefficient-parameter", "diverging-paths", "sparse", "generator"],
     variables = vars, parameters = String[], equations = eqs,
     generator = "eco(n) in $GENERATOR (factored form as in PHCpack Demo/eco5..eco8).",
     source = ["citation" => "A. Morgan, 'Solving polynomial systems using continuation for engineering and scientific problems', Prentice-Hall 1987, p. 148. PHCpack Demo/eco5..eco8 (this file: same formula at n = $n).",
               "url" => "https://github.com/janverschelde/PHCpack/blob/ed77cc0234aec45099d1d587a90a1617d5bfdacd/src/Demo/eco8",
               "file" => "PHCpack@ed77cc0234aec45099d1d587a90a1617d5bfdacd:src/Demo/eco8 (n = 8 member of the family); generated by $GENERATOR"],
     solutions = ["generic_nonsingular" => count,
                  "count_source" => "Mixed volume of the coefficient family = 2^(n-2) = $count (HC.jl 2.22.4 mixed_volume gives 2^(n-2) for every n = 8..16; PHCpack Demo/eco8: mixed volume 64 = 2^6, 64 solutions). See [verification] for the solve of the stored system."],
     parametrization = ["kind" => "coefficients", "description" => COEF_DESCRIPTION * " Generic count of the family = mixed volume = $count."],
     verification = verification)
end

# --------------------------------------------------------------------- noon

function noon(n)
    vars = ["x$i" for i in 1:n]
    eqs = [join(["x$i*x$j^2" for j in 1:n if j != i], " + ") * " - 1.1*x$i + 1" for i in 1:n]
    vars, eqs
end

function noon_system(n; verification)
    vars, eqs = noon(n)
    c = 3^n - 2n
    (name = "noon$n",
     description = "Noonburg's neural network (adaptive Lotka-Volterra), n = $n, the same family as tier_b/noon5 at a larger size: x_i * sum_{j != i} x_j^2 - 1.1 x_i + 1 = 0. $n cubics with $(n+1) terms each, full S_$n symmetry.",
     tags = ["large", "stress", "coefficient-parameter", "symmetric", "generator"],
     variables = vars, parameters = String[], equations = eqs,
     generator = "noon(n) in $GENERATOR (as PHCpack Demo/noon3..noon5, coefficient c = 1.1).",
     source = ["citation" => "V. W. Noonburg, 'A neural network modeled by an adaptive Lotka-Volterra system', SIAM J. Appl. Math. 49 (1989). PHCpack Demo/noon3..noon5 (this file: same formula at n = $n).",
               "url" => "https://github.com/janverschelde/PHCpack/blob/ed77cc0234aec45099d1d587a90a1617d5bfdacd/src/Demo/noon5",
               "file" => "PHCpack@ed77cc0234aec45099d1d587a90a1617d5bfdacd:src/Demo/noon5 (n = 5 member of the family); generated by $GENERATOR"],
     solutions = ["generic_nonsingular" => c,
                  "count_source" => "3^n - 2n = $c (PHCpack Demo/noon3..noon5: 21, 73, 233 = 3^n - 2n, equal to the mixed volume; HC.jl 2.22.4 mixed_volume gives 3^n - 2n for every n = 5..9). See [verification]."],
     parametrization = ["kind" => "coefficients", "description" => COEF_DESCRIPTION * " Generic count of the family = mixed volume = $c."],
     verification = verification)
end

# --------------------------------------------------------- discretized BVP

function bvp(N)
    h2 = (N + 1)^2
    vars = ["y$k" for k in 1:N]
    eqs = String[]
    for k in 1:N
        s = (k > 1 ? "y$(k-1) - " : "-") * "2*y$k" * (k < N ? " + y$(k+1)" : "")
        push!(eqs, s * " + (1/$h2)*lam*(1 + y$k^2) - (1/$h2)*f$k")
    end
    vars, ["lam"; ["f$k" for k in 1:N]], eqs
end

function bvp_system(N; verification)
    vars, params, eqs = bvp(N)
    (name = "bvp_quadratic_n$N",
     description = "Central finite-difference discretization of the two-point boundary value problem y'' = -lam (1 + y^2) + f(t), y(0) = y(1) = 0, on N = $N interior nodes t_k = k h, h = 1/$(N+1): y_{k-1} - 2 y_k + y_{k+1} + h^2 lam (1 + y_k^2) - h^2 f_k = 0 (y_0 = y_{N+1} = 0). Tridiagonal, 4-5 terms per equation, $N unknowns. With f = 0 this is the problem of Allgower-Bates-Sommese-Wampler (2006, Sec. 3.2); the continuous problem has two solutions for 0 < lam < lam*, one at lam*, none above. The forcing samples f_k are extra parameters so that the family has a provable generic root count. Intended use: tracking a few paths (e.g. the real solutions while lam or f moves) in a large sparse system; a full solve is impossible (2^$N solutions).",
     tags = ["large", "boundary-value", "sparse", "tridiagonal", "many-variables", "single-path", "parameter-homotopy", "generator"],
     variables = vars, parameters = params, equations = eqs,
     example_parameters = ["3"; fill("0", N)],
     generator = "bvp(N) in $GENERATOR: variables y1..yN, parameters lam, f1..fN, h^2 = 1/(N+1)^2 written as the literal 1/(N+1)^2.",
     source = ["citation" => "E. L. Allgower, D. J. Bates, A. J. Sommese, C. W. Wampler, 'Solution of polynomial systems derived from differential equations', Computing 76 (2006) 1-10, Sec. 3.2 (y'' = -lam(1 + y^2), y(0) = y(1) = 0); also W. Hao, J. D. Hauenstein, B. Hu, A. J. Sommese, 'A bootstrapping approach for computing multiple solutions of differential equations', J. Comput. Appl. Math. 258 (2014) 181-190. Discretization written for this collection (forcing term added).",
               "url" => "https://doi.org/10.1007/s00607-005-0132-4",
               "file" => "generated by $GENERATOR"],
     solutions = ["generic_nonsingular_formula" => "2^N = 2^$N",
                  "count_source" => "Bezout: the top-degree parts h^2 lam y_k^2 have no common nontrivial zero for lam != 0, so there are no solutions at infinity and exactly 2^N solutions counted with multiplicity for every f. For fixed lam != 0 the map y -> (residuals without the f terms) is finite onto C^N, so by generic smoothness all 2^N solutions are nonsingular for generic f. (2^$N does not fit in a TOML integer, hence the formula.) For f = 0 the count of nonsingular solutions is not claimed."],
     parametrization = ["kind" => "explicit", "description" => "lam: load parameter; f1..fN: forcing samples f(t_k). example_parameters = (lam = 3, f = 0), the unforced problem of Allgower et al."],
     verification = verification)
end

# ------------------------------------------------------- dense random systems

function dense(n, d)
    vars = ["x$i" for i in 1:n]
    ms = monomials(n, d)
    params = String[]
    eqs = String[]
    for i in 1:n
        terms = String[]
        for (k, e) in enumerate(ms)
            c = "c$(i)_$k"
            push!(params, c)
            push!(terms, term_string(c, monomial_string(vars, e)))
        end
        push!(eqs, join(terms, " + "))
    end
    vars, params, eqs
end

function dense_system(n, d)
    vars, params, eqs = dense(n, d)
    m = binomial(n + d, d)
    (name = "dense_n$(n)_d$(d)",
     description = "Dense square system: $n equations in $n unknowns, each a general polynomial of degree $d with all binomial($(n+d), $d) = $m monomials of degree <= $d, every coefficient an independent parameter ($(n*m) parameters, $(n*m) terms). Stresses evaluation and differentiation of many-term polynomials. Equation i is sum_k c{i}_{k} * m_k with the monomials m_k of degree <= $d in graded order (degree $d first, the constant term c{i}_$m last).",
     tags = ["large", "dense", "many-terms", "many-parameters", "total-degree", "coefficient-parameter", "generator"],
     variables = vars, parameters = params, equations = eqs,
     generator = "dense(n, d) in $GENERATOR: monomials of degree <= d in x1..xn ordered by total degree descending, then exponent vector lexicographically descending (x1^d first, constant last); parameter c{i}_{k} is the coefficient of the k-th monomial in equation i. A random instance: i.i.d. complex Gaussian coefficients.",
     source = ["citation" => "Bezout's theorem: a system of n generic polynomials of degree d in n unknowns has exactly d^n isolated nonsingular solutions (e.g. A. J. Sommese, C. W. Wampler, 'The Numerical Solution of Systems of Polynomials Arising in Engineering and Science', World Scientific 2005, Ch. 8, total-degree homotopy). Random dense systems of this type are the classical worst case for evaluation cost; system written for this collection.",
               "url" => "https://doi.org/10.1142/5763",
               "file" => "generated by $GENERATOR"],
     solutions = ["generic_nonsingular" => d^n,
                  "count_source" => "Bezout number d^n = $(d)^$(n) = $(d^n), attained for generic coefficients (the family contains x_i^$d - 1 = 0, which has d^n nonsingular solutions and none at infinity)."],
     parametrization = ["kind" => "explicit", "description" => "Every monomial coefficient is a parameter (this is the coefficient-parameter family itself). Total-degree start system x_i^$d - 1 is the member with c{i}_{k} = 1 for the pure power x_i^$d, -1 for the constant, 0 otherwise."],
     verification = ["status" => "not-needed", "method" => "count is the Bezout number of a generic dense system"])
end

# ------------------------------------------------- power flow, complete graph

function powerflow(n)  # n + 1 buses 0..n, bus 0 is the slack bus
    vars = String[]
    for i in 1:n
        push!(vars, "x$i", "y$i")
    end
    bname(i, j) = "b$(min(i, j))$(max(i, j))"
    params = [[bname(i, j) for i in 0:n for j in i+1:n]; ["P$i" for i in 1:n]]
    eqs = String[]
    for i in 1:n
        terms = ["b0$i*y$i"; ["$(bname(i, j))*(y$i*x$j - x$i*y$j)" for j in 1:n if j != i]]
        push!(eqs, join(terms, " + ") * " - P$i")
    end
    append!(eqs, ["x$i^2 + y$i^2 - 1" for i in 1:n])
    vars, params, eqs
end

function powerflow_system(n; verification)
    @assert n <= 9 # single-digit bus indices in the parameter names
    vars, params, eqs = powerflow(n)
    c = binomial(2n, n)
    (name = "lossless_powerflow_k$(n+1)",
     description = "Lossless PV-bus power flow (equivalently Kuramoto equilibria) on the complete graph K$(n+1), the same formulation as tier_a/lossless_powerflow_k4 at a larger size. Bus 0 is the slack/reference bus with angle 0 (x0 = 1, y0 = 0); for buses i = 1..$n, x_i = cos(theta_i), y_i = sin(theta_i). Active power balance: sum_j b_ij sin(theta_i - theta_j) = P_i, with sin(theta_i - theta_j) = y_i x_j - x_i y_j. Parameters: $(binomial(n+1, 2)) line susceptances b_ij and $n injections P_i.",
     tags = ["large", "power-flow", "kuramoto", "parameter-homotopy", "generator"],
     variables = vars, parameters = params, equations = eqs,
     generator = "powerflow(n) in $GENERATOR for n + 1 buses (n <= 9 because the parameter names b_ij use single-digit bus indices).",
     source = ["citation" => "J. Baillieul, C. I. Byrnes, 'Geometric critical point analysis of lossless power system models', IEEE Trans. Circuits Syst. 29(11), 1982 (bound C(2n-2, n-1) for n buses); T. Chen, R. Davis, D. Mehta, 'Counting equilibria of the Kuramoto model using birationally invariant intersection index', SIAM J. Appl. Algebra Geom. 2(4), 2018. Formulation written for this collection (as tier_a/lossless_powerflow_k4).",
               "url" => "https://doi.org/10.1109/TCS.1982.1085093",
               "file" => "generated by $GENERATOR"],
     solutions = ["generic_nonsingular" => c,
                  "count_source" => "Baillieul-Byrnes bound C(2N-2, N-1) = C($(2n), $n) = $c for N = $(n+1) buses, attained for generic parameters on the complete graph (Chen-Davis-Mehta 2018). See [verification]."],
     parametrization = ["kind" => "explicit"],
     verification = verification)
end

# ------------------------------------------- lines on a septic fivefold

# Coefficients of t^k of F(a1 t + b1, ..., a4 t + b4, t) for a dense F of
# degree 7 in x1..x5 with parameter coefficients q1..q791 (constant term 1).
function fano_septic()
    nx, d = 5, 7
    ms = monomials(nx, d)        # constant monomial is last
    vars = [["a$i" for i in 1:nx-1]; ["b$i" for i in 1:nx-1]]
    params = ["q$k" for k in 1:length(ms)-1]
    # eqs[k+1]: Dict(exponent vector over vars => Dict(q index (0 = constant 1) => integer coefficient))
    eqs = [Dict{Vector{Int},Dict{Int,Int}}() for _ in 0:d]
    for (qi, e) in enumerate(ms)
        qidx = qi == length(ms) ? 0 : qi
        for ks in Iterators.product((0:e[i] for i in 1:nx-1)...)
            tdeg = sum(ks) + e[nx]
            mono = [collect(ks); [e[i] - ks[i] for i in 1:nx-1]]
            coef = prod(binomial(e[i], ks[i]) for i in 1:nx-1)
            dct = get!(eqs[tdeg+1], mono, Dict{Int,Int}())
            dct[qidx] = get(dct, qidx, 0) + coef
        end
    end
    strs = String[]
    for k in 0:d
        terms = String[]
        for mono in sort(collect(keys(eqs[k+1])); by = m -> (sum(m), reverse(m)))
            for (qidx, coef) in sort(collect(eqs[k+1][mono]))
                c = qidx == 0 ? "" : "q$qidx"
                m = monomial_string(vars, mono)
                t = term_string(c, m)
                t == "" && (t = "1")
                push!(terms, coef == 1 ? t : "$coef*$t")
            end
        end
        push!(strs, join(terms, " + "))
    end
    vars, params, strs
end

function fano_septic_system(; verification)
    vars, params, eqs = fano_septic()
    (name = "fano_septic_lines",
     description = "Lines on a generic septic (degree 7) hypersurface in P^5, the next member after hcjl/fano_quintic_lines (lines on a quintic threefold, 2875). F is a dense affine septic in x1..x5 with all binomial(12, 7) = 792 monomial coefficients; the constant coefficient is set to 1, leaving 791 parameters q1..q791. A line is x(t) = (a1 t + b1, a2 t + b2, a3 t + b3, a4 t + b4, t); the 8 equations are the coefficients of t^0..t^7 of F(x(t)) (equation k+1 is the coefficient of t^k), in 8 unknowns a1..a4, b1..b4. 11440 terms in total.",
     tags = ["large", "enumerative-geometry", "parameter-homotopy", "many-terms", "many-parameters", "large-count", "generator"],
     variables = vars, parameters = params, equations = eqs,
     generator = "fano_septic() in $GENERATOR: monomials of F ordered by total degree descending, then exponent vector (x1..x5) lexicographically descending; q_k is the coefficient of the k-th monomial, the constant monomial (last) has coefficient 1. Substituting x_i = a_i t + b_i (i <= 4), x5 = t and expanding with the binomial theorem gives the stored equations (checked numerically against direct evaluation of F(x(t)) at random points).",
     source = ["citation" => "Number of lines on the generic hypersurface of degree 2n-1 in P^(n+1): 1, 27, 2875, 698005, ... (OEIS A027363; B. L. van der Waerden, 'Zur algebraischen Geometrie' series). Formulation as in HC.jl fano_quintic() (test/test_systems.jl:730-870), extended to degree 7 for this collection.",
               "url" => "https://oeis.org/A027363",
               "file" => "generated by $GENERATOR; formulation of HomotopyContinuation.jl@0cbf1e27d062d7eb01d962c3b948c39444c9a62c:test/test_systems.jl:730-870 (fano_quintic)"],
     solutions = ["generic_nonsingular" => 698005,
                  "count_source" => "Schubert-calculus count: a generic septic hypersurface in P^5 contains exactly 698005 lines (OEIS A027363, term n = 4), all reduced. For generic coordinates all of them lie in the chart x(t) = (a t + b, t)."],
     parametrization = ["kind" => "explicit"],
     verification = verification)
end

# ------------------------------------------------------------------ systems

const HC = "2.22.4"
const DATE = "2026-11-08"

function main()
    write_system(katsura_system(15; verification = [
        "status" => "verified",
        "method" => "HC.jl solve (polyhedral, HC's affine variant) of the stored system: 32768 paths, 32768 nonsingular (14 s, 16 threads). HC.jl mixed_volume (torus) = 32730. The generator reproduces tier_b/katsura12 exactly for n = 12.",
        "hc_version" => HC, "verified_count" => 32768, "date" => DATE]))
    write_system(katsura_system(20; verification = [
        "status" => "not-verified",
        "method" => "not solved here (1048576 paths). Same generator as katsura15 (verified, 32768) and tier_b/katsura12."]))
    write_system(cyclic_system(10; mv = 35940, verification = [
        "status" => "verified",
        "method" => "HC.jl mixed_volume = 35940; HC.jl solve (polyhedral) with independent complex Gaussian coefficients on every monomial: 35940 paths, 35940 nonsingular (7 s, 16 threads). Equations equal PHCpack Demo/cyclic10 (difference 0 at a random point).",
        "hc_version" => HC, "verified_count" => 35940, "date" => DATE]))
    write_system(eco_system(15; count = 8192, verification = [
        "status" => "verified",
        "method" => "HC.jl solve (polyhedral) of the stored system: 8192 paths, 8192 nonsingular (16 threads). HC.jl mixed_volume of the support = 8192. With i.i.d. complex Gaussian coefficients on every monomial the same solve returns 7919 nonsingular + 271 singular (the random member is badly conditioned), so the family count rests on the mixed volume. The generator reproduces tier_b/eco8 exactly for n = 8.",
        "hc_version" => HC, "verified_count" => 8192, "date" => DATE]))
    write_system(noon_system(9; verification = [
        "status" => "verified",
        "method" => "HC.jl solve (polyhedral): stored system 19665 paths, 19665 nonsingular; independent complex Gaussian coefficients on every monomial: 19665 nonsingular. HC.jl mixed_volume = 19665. The generator reproduces tier_b/noon5 exactly for n = 5.",
        "hc_version" => HC, "verified_count" => 19665, "date" => DATE]))
    write_system(bvp_system(100; verification = [
        "status" => "partial",
        "method" => "Count by the argument in count_source. Small members checked with HC.jl total-degree solve at lam = 3: N = 8 gives 256 nonsingular for random complex f and for f = 0; N = 10 gives 1024 nonsingular for both.",
        "hc_version" => HC, "date" => DATE]))
    write_system(dense_system(6, 4))
    write_system(dense_system(8, 4))
    write_system(dense_system(4, 10))
    write_system(powerflow_system(9; verification = [
        "status" => "verified",
        "method" => "HC.jl monodromy_solve from a random start pair found by HC (find_start_pair), max_loops_no_progress = 10: 48620 solutions (204 s, 16 threads). HC.jl mixed_volume at random complex parameters = 257536. The generator reproduces tier_a/lossless_powerflow_k4 exactly for n = 3.",
        "hc_version" => HC, "verified_count" => 48620, "date" => DATE]))
    write_system(fano_septic_system(; verification = [
        "status" => "not-verified",
        "method" => "Not solved (698005 paths). The expanded equations were checked against direct evaluation of F(x(t)) = sum_k eq_k t^k at random complex a, b, q, t (relative error 1e-14)."]))
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
