using Vaje08
using Test

@testset "Baza" begin
  @test isapprox(1, Vaje08.h00(0))
  @test isapprox(0, Vaje08.h00(1))
  @test isapprox(0, Vaje08.h01(0))
  @test isapprox(1, Vaje08.h01(1))
  @test isapprox(0, Vaje08.h10(0))
  @test isapprox(0, Vaje08.h10(1))
  @test isapprox(0, Vaje08.h11(0))
  @test isapprox(0, Vaje08.h11(1))
end

@testset "Hermitova interpolacija" begin
  # Funkcija: x^3
  xint = [1, 3]
  y = xint .^ 3
  dy = 3 * xint .^ 2

  @test isapprox(8, hermiteInt(2, xint, y, dy))
end

@testset "Hermitov zlepek" begin
  # Funkcija: x^3
  x = [0, 2, 3, 5]
  y = x .^ 3
  dy = 3 * x .^ 2
  Z = HermitovZlepek(x, y, dy)

  @test isapprox(1, vrednost(1, Z))
  @test isapprox(0, vrednost(0, Z))
  @test isapprox(2.5^3, vrednost(2.5, Z))
  @test_throws BoundsError vrednost(-1, Z)
end