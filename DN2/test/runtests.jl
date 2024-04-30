using DN2
using Test


@testset "normal_CDF_test" begin
  tol = 10^-10

  @test abs(DN2.normal_CDF(0.0) - 0.5) <= tol
  @test abs(DN2.normal_CDF(0.5) - 0.6914624612925835) <= tol
  @test abs(DN2.normal_CDF(-0.5) - 0.3085375387074165) <= tol
  @test abs(DN2.normal_CDF(1.0) - 0.8413447460843065) <= tol

  @test abs(DN2.normal_CDF(10) - 1) <= tol
  @test abs(DN2.normal_CDF(-10)) <= tol
end


@testset "bin_coef" begin
  @test DN2.binomial_coefficient(3, 2) == 3
  @test DN2.binomial_coefficient(3, 1) == 3
  @test DN2.binomial_coefficient(10, 5) == 252

end

@testset "bezier_linear" begin
  control_points = [(0, 0), (1, 1)]
  @test DN2.bezier(control_points, 0.0) == (0, 0)
  @test DN2.bezier(control_points, 0.15) == (0.15, 0.15)
  @test DN2.bezier(control_points, 0.5) == (0.5, 0.5)
  @test DN2.bezier(control_points, 0.75) == (0.75, 0.75)
  @test DN2.bezier(control_points, 1.0) == (1.0, 1.0)
end

@testset "bezier_quadratic" begin
  control_points = [(0, 0), (1, 1), (2, 0)]
  @test DN2.bezier(control_points, 0.0) == (0, 0)
  @test DN2.bezier(control_points, 0.5) == (1.0, 0.5)
  @test DN2.bezier(control_points, 1.0) == (2.0, 0.0)
end
