mutable struct Statistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    nparametergenerations::Int
    nsolutions_development::Vector{Int}
end
Statistics(nsolutions::Int) = Statistics(0, 0, 1, [nsolutions])

function trackedpath!(stats::Statistics, retcode)
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

function finished!(stats, nsolutions)
    push!(stats.nsolutions_development, nsolutions)
end
