using HomotopyContinuation
using TestSystems

F = equations(katsura8())

@time solve(F)

F
