using DN2
using Test
using FastGaussQuadrature


function accurate_gauss_CDF(x)
  f(x) = exp(-x^2 / 2) / sqrt(2 * pi)
  a = -10
  b = x
  n = 60
  xs, ws = gausslegendre(n)
  xs = (b - a) / 2 * xs .+ (b + a) / 2
  ws = (b - a) / 2 * ws
  return sum(ws .* f.(xs))
end


@testset "gaussian_CDF_test" begin
  tol = 10^-10

  @test abs(DN2.gaussian_CDF(0.0) - 0.5) / 0.5 <= tol
  @test abs(DN2.gaussian_CDF(0.5) - 0.6914624612925835) / 0.6914624612925835 <= tol
  @test abs(DN2.gaussian_CDF(-0.5) - 0.3085375387074165) / 0.3085375387074165 <= tol
  @test abs(DN2.gaussian_CDF(1.0) - 0.8413447460843065) / 0.8413447460843065 <= tol

  @test abs(DN2.gaussian_CDF(10) - 1) / 1 <= tol
  @test abs(DN2.gaussian_CDF(-10)) <= tol
  @test abs(DN2.gaussian_CDF(-100000) - 0) <= tol
end


@testset "gaussian_CDF_test_range" begin
  tol = 10^-10

  for val in range(-1, stop=1, length=200)
    true_val = accurate_gauss_CDF(val)
    @test abs(DN2.gaussian_CDF(val) - true_val) / abs(true_val) <= tol
  end
end


@testset "binomial_coefficient" begin
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
