using DN1
using Test

@testset "nakljucno izbrane tocke" begin
  x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  y = [1.0, 3.0, 1.0, 2.0, 0.0, 6.0]
  s = interpoliraj(x, y)
  p1 = vrednost(s, 5)
  p2 = vrednost(s, 4.9999)
  p3 = vrednost(s, 0.0)
  p4 = vrednost(s, 0.5)

  @test p1 == 6.0
  @test abs(p2 - 5.99) < 0.01
  @test p3 == 1.0
  @test p4 >= 1.0 && p4 <= 3.0
end

@testset "sinus" begin
  x = range(0, stop=2 * pi, length=5)
  y = sin.(x)
  s = interpoliraj(x, y)

  p1 = vrednost(s, 0)
  p2 = vrednost(s, pi / 4)
  p3 = vrednost(s, pi / 2)
  p4 = vrednost(s, pi)
  p5 = vrednost(s, 2 * pi)

  @test p1 == 0.0
  @test p2 - sqrt(2) / 2 < 1e-15
  @test p3 ≈ 1.0
  @test abs(p4) < 1e-15
  @test abs(p5) < 1e-15
end