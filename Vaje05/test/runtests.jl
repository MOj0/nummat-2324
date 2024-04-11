using Vaje05
using Test

@testset "3x3 prehodna matrika" begin
  P = [0.1 0.4 0.5; 0 0.5 0.5; 0.2 0 0.8]
  v, lambda = Vaje05.potencna(P, [-1, 0, 0], 100, 1e-10)

  @test lambda ≈ 1
  @test v ≈ [1, 1, 1]
end