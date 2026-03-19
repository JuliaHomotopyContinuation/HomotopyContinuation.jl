using JET
using Test
using HomotopyContinuation

# Filter known reports: cos/sin/sqrt(::IComplexF64) are not implemented because
# adding them to Base causes +4s TTFX regression via @generated execute_instructions! recompilation.
# These only trigger when certifying systems with transcendental functions.
const KNOWN_MISSING = Set(["cos", "sin", "sqrt"])

function is_known_icomplexf64_report(r)
    s = sprint(show, r)
    occursin("IComplexF64", s) &&
        any(f -> occursin("no matching method found `$f(", s), KNOWN_MISSING)
end

@testset "JET checks" begin
    rep = JET.report_package(HomotopyContinuation; target_modules = (HomotopyContinuation,))
    reports = JET.get_reports(rep)
    filtered = filter(!is_known_icomplexf64_report, reports)
    @show rep
    if length(filtered) != length(reports)
        @info "Filtered $(length(reports) - length(filtered)) known IComplexF64 reports"
    end
    @test length(filtered) == 0
end
