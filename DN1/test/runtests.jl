using DN1
using Test

@testset "koeficienti_2_tocki" begin
  x = [0.0, 1.0]
  y = [1.0, 10.0]
  s = interpoliraj(x, y)

  @test s.coefficients[1] ≈ 1.0
  @test s.coefficients[2] ≈ 9.0
  @test s.coefficients[3] ≈ 0
  @test s.coefficients[4] ≈ 0

  @test vrednost(s, 0.5) ≈ (10 - 1) / 2 + 1
  @test vrednost(s, 1 / 3) ≈ (10 - 1) / 3 + 1
end

@testset "koeficienti_3_tocke" begin
  x = [0.0, 1.0, 2.0]
  y = [1.0, 5.0, 2.0]
  s = interpoliraj(x, y)

  @test s.coefficients ≈ [1.0 5.75 0.0 -1.75; 5.0 0.5 -5.25 1.75]

  @test vrednost(s, 0.5) ≈ 3.65625
  @test vrednost(s, 0.75) ≈ 4.57421875
  @test vrednost(s, 0.95) ≈ 4.96209375
  @test vrednost(s, 1.25) ≈ 4.82421875
  @test vrednost(s, 1.7325) ≈ 3.2371277304687505
end

@testset "Zlepek na 6 tockah" begin
  x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  y = [1.0, 3.0, 1.0, 2.0, 0.0, 6.0]
  s = interpoliraj(x, y)
  p1 = vrednost(s, 5)
  p2 = vrednost(s, 4.9999)
  p3 = vrednost(s, 0.0)
  p4 = vrednost(s, 0.5)

  @test p1 ≈ 6.0
  @test abs(p2 - 5.99) < 0.01
  @test p3 ≈ 1.0
  @test p4 >= 1.0 && p4 <= 3.0
end

@testset "sinus" begin
  x = range(0, stop=2 * pi, length=50)
  y = sin.(x)
  s = interpoliraj(x, y)

  diff = [abs(vrednost(s, xi) - sin(xi)) for xi in x]
  for d in diff
    @test d < 1e-15
  end
end