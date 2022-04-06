using BenchmarkTools
using PrettyTables

function judge_results(result1_name, result2_name)
    r1 = BenchmarkTools.load(result1_name)[1]
    r2 = BenchmarkTools.load(result2_name)[1]
    res = BenchmarkTools.judge(median(r1), median(r2))

    println("Interpreted:")

    data = vcat(
        map(res["InterpretedSystem"]) do (system, group)
            [system group["evaluate"] group["evaluate_and_jacobian"] group["taylor_1"] group["taylor_2"]]
        end...,
    )

    pretty_table(
        data;
        header = ["System", "evaluate", "evaluate_and_jacobian", "taylor_1", "taylor_2"],
    )
end

