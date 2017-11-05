@testset "Patch" begin
  A = big.(rand(4, 4))
  b = big.(rand(4))
  c = copy(b)

  LU = lufact(A)
  HomotopyContinuation.my_A_ldiv_B!(LU, b)

  @test b ≈ (A \ c)
end
