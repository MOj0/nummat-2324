using Test, Vaje09

@testset "Trapezno pravilo" begin
  x, w = trapez(1, 4, 3)
  @test isapprox(x, [1, 2, 3, 4])
  @test isapprox(w, [0.5, 1, 1, 0.5])
end

@testset "Simpsonovo pravilo" begin
  x, w = simpson(1, 7, 3)
  @test isapprox(x, [1, 2, 3, 4, 5, 6, 7])
  @test isapprox(w, 1 / 3 * [1, 4, 2, 4, 2, 4, 1])
end