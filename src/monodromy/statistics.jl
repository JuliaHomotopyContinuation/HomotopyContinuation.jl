mutable struct Statistics
    ntrackedpaths::Int
    ntrackingfailures::Int
    nreal::Int
    nparametergenerations::Int
    nsolutions_development::Vector{Int}
end

Statistics(nsolutions::Int) = Statistics(0, 0, 0, 1, [nsolutions])
function Statistics(solutions)
    stats = Statistics(length(solutions))
    for s in solutions
        checkreal!(stats, s)
    end
    stats
end


Base.show(io::IO, S::Statistics) = Utilities.print_fieldnames(io, S)
Base.show(io::IO, ::MIME"application/prs.juno.inline", S::Statistics) = S


function trackedpath!(stats::Statistics, retcode)
    if retcode == :success
        stats.ntrackedpaths += 1
    else
        stats.ntrackingfailures += 1
    end
end


function checkreal!(stats::Statistics, y)
    if Utilities.isrealvector(y)
        stats.nreal +=1
    end
end

function generated_parameters!(stats::Statistics, nsolutions::Int)
    stats.nparametergenerations += 1
    push!(stats.nsolutions_development, nsolutions)
end

function finished!(stats, nsolutions)
    push!(stats.nsolutions_development, nsolutions)
end
