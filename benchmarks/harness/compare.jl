# Compare the endpoints of harness runs path by path with a reference run (ideally one
# whose endpoints were certified, e.g. a run of the commit the change starts from).
#   julia --project=benchmarks/harness benchmarks/harness/compare.jl <sha> <label-prefix> <reference> [out.tsv]
# Every run `<sha>_<label-prefix>*_<mode>[_scaled]` is compared with the reference run of
# the same mode (`<reference>_<mode>[_scaled]`, e.g. reference `62742441_baseline`).
# Per path i whose reference run succeeded:
#   ok       reldist(x, gt[i]) < 1e-6
#   visible  wrong and the endpoint coincides with another successful endpoint of this run
#   silent   wrong, equals some other reference endpoint, and no duplicate in this run
#   unknown  wrong and matches no reference endpoint (new point, or garbage)
#   fail     no success (paths where the reference failed are skipped)
using Serialization, Printf
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EP = joinpath(ROOT, "benchmarks", "cache", "endpoints")
reldist(a, b) = maximum(abs, a .- b) / max(1.0, maximum(abs, b))

# indices j with reldist(x, xs[j]) < tol, using a sort on real(x[1])
struct Index; xs; perm; keys; end
Index(xs) = (p = sortperm(xs; by = x -> real(x[1])); Index(xs, p, [real(xs[i][1]) for i in p]))
function matches(ix::Index, x; tol = 1e-6)
    w = 2tol * max(1.0, maximum(abs, x))
    lo = searchsortedfirst(ix.keys, real(x[1]) - w); hi = searchsortedlast(ix.keys, real(x[1]) + w)
    [ix.perm[j] for j in lo:hi if reldist(x, ix.xs[ix.perm[j]]) < tol]
end

length(ARGS) ≥ 3 || error("usage: compare.jl <sha> <label-prefix> <reference> [out.tsv]")
sha, prefix, reference = ARGS[1], ARGS[2], ARGS[3]
out = open(length(ARGS) ≥ 4 ? ARGS[4] : "compare_$(sha)_$(prefix).tsv", "w")
println(out, join(["variant", "scaled", "system", "paths", "gt_ok", "ok", "fail", "visible", "silent", "unknown"], '\t'))
for f in sort(readdir(EP))
    m = match(r"^([\w-]+?)_([^_]+)_(interpreted|compiled)(_scaled)?\.jls$", f)
    m === nothing && continue
    m.captures[1] == sha && startswith(m.captures[2], prefix) || continue
    variant = m.captures[2]; mode = m.captures[3]; sc = m.captures[4] !== nothing
    ref = joinpath(EP, "$(reference)_$(mode)$(sc ? "_scaled" : "").jls")
    isfile(ref) || (println("skipped $f: no reference run $(basename(ref))"); continue)
    run = deserialize(joinpath(EP, f))
    gt = deserialize(ref)
    for name in sort(collect(keys(run)))
        g = gt[name]; r = run[name]
        gsucc = findall(==(:success), g.status)
        gix = Index(g.x[gsucc])
        rsucc = findall(==(:success), r.status)
        rix = Index(r.x[rsucc])
        ok = fail = vis = sil = unk = 0
        for i in gsucc
            if r.status[i] !== :success; fail += 1; continue; end
            x = r.x[i]
            if reldist(x, g.x[i]) < 1e-6; ok += 1; continue; end
            if length(matches(rix, x)) > 1; vis += 1
            elseif isempty(matches(gix, x)); unk += 1
            else sil += 1; end
        end
        println(out, join(Any[variant, sc, name, length(r.status), length(gsucc), ok, fail, vis, sil, unk], '\t'))
    end
end
close(out)
