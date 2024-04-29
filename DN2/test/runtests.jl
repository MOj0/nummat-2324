using DN2
using Test


@testset "normal_CDF_test" begin
  tol = 10^-10

  @test abs(normal_CDF(0.5) - 0.6914624612925835) <= tol
  @test abs(normal_CDF(0.0) - 0.5) <= tol
  @test abs(normal_CDF(1.0) - 0.8413447460843065) <= tol
end
