@testset "Compiled cache thread safety" begin
    if Threads.nthreads() < 2
        @test_skip "requires at least two Julia threads"
    else
        n_compiled_cache_test_systems = 128

        @var x t
        systems = [
            System([x^2 + (10_000 + i) * x + i], [x]) for
            i = 1:n_compiled_cache_test_systems
        ]
        homotopies = [
            Homotopy([x^2 + (20_000 + i) * x + t + i], [x], t) for
            i = 1:n_compiled_cache_test_systems
        ]

        function compile_concurrently(inputs, compile)
            ready = Threads.Atomic{Int}(0)
            start = Base.Event()
            tasks = map(1:2) do worker
                Threads.@spawn begin
                    Threads.atomic_add!(ready, 1)
                    wait(start)
                    compiled = Any[]
                    for i = worker:2:length(inputs)
                        c = compile(inputs[i])
                        push!(compiled, (i, c, ModelKit.interpret(typeof(c))))
                    end
                    compiled
                end
            end
            while ready[] < length(tasks)
                yield()
            end
            notify(start)
            reduce(vcat, fetch.(tasks))
        end

        function run_while_lock_is_held(f, cache_lock)
            task = lock(cache_lock) do
                task = Threads.@spawn f()
                @test Base.timedwait(() -> istaskdone(task), 0.25) == :timed_out
                task
            end
            fetch(task)
        end

        # Warm up both constructor paths before checking that they wait for the
        # corresponding cache lock.
        CompiledSystem(System([x^2 + 30_001 * x + 1], [x]))
        CompiledHomotopy(Homotopy([x^2 + 40_001 * x + t + 1], [x], t))

        blocking_system = System([x^2 + 30_002 * x + 2], [x])
        compiled_system = run_while_lock_is_held(
            () -> CompiledSystem(blocking_system),
            ModelKit._TSYSTEM_TABLE_LOCK,
        )
        cached_system = run_while_lock_is_held(
            () -> ModelKit.interpret(typeof(compiled_system)),
            ModelKit._TSYSTEM_TABLE_LOCK,
        )
        @test cached_system([0.25]) ≈ blocking_system([0.25])

        blocking_homotopy = Homotopy([x^2 + 40_002 * x + t + 2], [x], t)
        compiled_homotopy = run_while_lock_is_held(
            () -> CompiledHomotopy(blocking_homotopy),
            ModelKit._THOMOTOPY_TABLE_LOCK,
        )
        cached_homotopy = run_while_lock_is_held(
            () -> ModelKit.interpret(typeof(compiled_homotopy)),
            ModelKit._THOMOTOPY_TABLE_LOCK,
        )
        @test cached_homotopy([0.25], 0.5) ≈ blocking_homotopy([0.25], 0.5)

        system_entries_before = sum(length, values(ModelKit.TSYSTEM_TABLE); init = 0)
        compiled_systems = compile_concurrently(systems, CompiledSystem)
        system_entries_after = sum(length, values(ModelKit.TSYSTEM_TABLE); init = 0)

        @test length(compiled_systems) == n_compiled_cache_test_systems
        @test system_entries_after - system_entries_before == n_compiled_cache_test_systems
        @test all(compiled_systems) do (i, compiled, cached)
            compiled([0.25]) ≈ cached([0.25]) ≈ systems[i]([0.25])
        end
        @test allunique(typeof(c) for (_, c, _) in compiled_systems)

        homotopy_entries_before = sum(length, values(ModelKit.THOMOTOPY_TABLE); init = 0)
        compiled_homotopies = compile_concurrently(homotopies, CompiledHomotopy)
        homotopy_entries_after = sum(length, values(ModelKit.THOMOTOPY_TABLE); init = 0)

        @test length(compiled_homotopies) == n_compiled_cache_test_systems
        @test homotopy_entries_after - homotopy_entries_before ==
              n_compiled_cache_test_systems
        @test all(compiled_homotopies) do (i, compiled, cached)
            compiled([0.25], 0.5) ≈ cached([0.25], 0.5) ≈ homotopies[i]([0.25], 0.5)
        end
        @test allunique(typeof(c) for (_, c, _) in compiled_homotopies)
    end
end
