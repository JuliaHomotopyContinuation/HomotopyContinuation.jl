@inline function tmap(f::F, solver::S, src; blocksize=20) where {F<:Function, S}
    n = length(src)
    dst₁ = f(solver, 1, src[1])
    dst = Vector{typeof(dst₁)}(undef, n)
    dst[1] = dst₁
    for k=2:n
        dst[k] = f(solver, 1, src[k])
    end
    dst
end

@inline function tmap(f::F, solvers::Vector{S}, src; blocksize=20) where {F<:Function, S}
    n = length(src)
    dst₁ = f(solvers[1], 1, src[1])
    dst = Vector{typeof(dst₁)}(undef, n)
    dst[1] = dst₁
    i = Threads.Atomic{Int}(2)
    blocksize = min(blocksize, div(blocksize, length(solvers)))

    function threads_fun()
        tid = Threads.threadid()
        solver = solvers[tid]
        thread_loop!(i, n, blocksize) do k
            dst[k] = f(solver, tid, src[k])
        end
    end

    ccall(:jl_threading_run, Ref{Cvoid}, (Any,), threads_fun)

    dst
end

function tforeach(f::F, solver, src; blocksize=20) where {F<:Function}
    foreach(x -> f(solver, 1, x), src)
end

function tforeach(f::F, solvers::AbstractVector, src; blocksize=20) where {F<:Function}
    i = Threads.Atomic{Int}(1)
    blocksize = min(blocksize, div(blocksize, length(solvers)))
    function threads_fun()
        tid = Threads.threadid()
        solver = solvers[tid]
        n = length(src)
        thread_loop!(i, n, blocksize) do k
            f(solver, tid, src[k])
        end
    end

    ccall(:jl_threading_run, Ref{Cvoid}, (Any,), threads_fun)

    nothing
end

function thread_loop!(f::F, i, n, blocksize) where {F<:Function}
    while true
        k = Threads.atomic_add!(i, blocksize)
        if k > n
            break
        end
        upper = min(k + blocksize, n)
        while k ≤ upper
            f(k)
            k += 1
        end
    end
end
