mutable struct Statistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    nparametergenerations::Int
    nsolutions_development::Vector{Int}
end
Statistics() = Statistics(0, 0, 0, Int[])

function pathtracked!(stats::Statistics, retcode)
    if retcode == :success
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end

function generated_parameters!(stats::Statistics, nsolutions::Int)
    stats.nparametergenerations += 1
    push!(stats.nsolutions_development, nsolutions)
end
