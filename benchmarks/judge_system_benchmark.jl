using BenchmarkTools
using PrettyTables

function judge_results(result1_name, result2_name)
    r1 = minimum(BenchmarkTools.load(result1_name)[1])
    r2 = minimum(BenchmarkTools.load(result2_name)[1])
    res = BenchmarkTools.judge(r1, r2)

    println("Interpreted [PR / main]:")
    data = vcat(
        map(res["InterpretedSystem"]) do (system, group)
            [system group["evaluate"].ratio.time group["evaluate_and_jacobian"].ratio.time group["taylor_1"].ratio.time group["taylor_2"].ratio.time]
        end...,
    )

    pretty_table(
        data;
        header = ["System", "evaluate", "evaluate_and_jacobian", "taylor_1", "taylor_2"],
    )


    println("Compiled [PR / main]:")
    data = vcat(
        map(res["CompiledSystem"]) do (system, group)
            [system group["evaluate"].ratio.time group["evaluate_and_jacobian"].ratio.time]
        end...,
    )
    pretty_table(data; header = ["System", "evaluate", "evaluate_and_jacobian"])


    println("Compiled / Interpreted [main]:")
    j = judge(r2["InterpretedSystem"], r2["CompiledSystem"])
    data = vcat(
        map(j) do (system, group)
            [system group["evaluate"].ratio.time group["evaluate_and_jacobian"].ratio.time]
        end...,
    )
    pretty_table(data; header = ["System", "evaluate", "evaluate_and_jacobian"])

    println("Compiled / Interpreted [PR]:")
    j = judge(r1["InterpretedSystem"], r1["CompiledSystem"])
    data = vcat(
        map(j) do (system, group)
            [system group["evaluate"].ratio.time group["evaluate_and_jacobian"].ratio.time]
        end...,
    )
    pretty_table(data; header = ["System", "evaluate", "evaluate_and_jacobian"])
    return r1, r2
end
