mutable struct MonodromyStatistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    nreal::Int
    nparametergenerations::Int
    nsolutions_development::Vector{Int}
end

MonodromyStatistics(nsolutions::Int) = Statistics(0, 0, 0, 1, [nsolutions])
function MonodromyStatistics(solutions)
    stats = MonodromyStatistics(length(solutions))
    for s in solutions
        checkreal!(stats, s)
    end
    stats
end


Base.show(io::IO, S::MonodromyStatistics) = print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::MonodromyStatistics) = S


function trackedpath!(stats::MonodromyStatistics, retcode)
    if retcode == PathTrackerStatus.success
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end


function checkreal!(stats::MonodromyStatistics, y)
    if isrealvector(y)
        stats.nreal +=1
    end
end

function generated_parameters!(stats::MonodromyStatistics, nsolutions::Int)
    stats.nparametergenerations += 1
    push!(stats.nsolutions_development, nsolutions)
end

function finished!(stats, nsolutions)
    push!(stats.nsolutions_development, nsolutions)
end
