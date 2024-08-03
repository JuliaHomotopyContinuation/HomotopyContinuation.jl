using SnoopCompile, HomotopyContinuation

inf_timing = @snoopi include(joinpath(@__DIR__, "snoop.jl"))
pc = SnoopCompile.parcel(inf_timing)

dir = mktempdir()
@info "Writing precompile results into $dir"

SnoopCompile.write(dir, pc)

mv(
    joinpath(dir, "precompile_HomotopyContinuation.jl"),
    joinpath(@__DIR__, "../src/precompile.jl");
    force = true,
)
