using Vaje07
using Test

@testset "Presek dveh kroznic" begin
  # Presek kroznic (x-1)^2 + y^2 = 4 in (x+1)^2 + y^2 = 4
  F(x) = [(x[1] - 1)^2 + x[2]^2 - 4, (x[1] + 1)^2 + x[2]^2 - 4]

  JF(x) = 2 * [(x[1]+1) x[2]; x[1]-1 x[2]]
  fdf(x) = (F(x), JF(x))
  x, it = newton(fdf, [0, 1])
  @test isapprox(F(x), zeros(2))
  @test isapprox(x, [0, sqrt(3)])
end