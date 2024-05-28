using DN3
using Test
using DifferentialEquations


function exact_sol(l, t, theta0, dtheta0)
  function pendulum!(du, u, p, t)
    g, l = p
    du[1] = u[2]
    du[2] = -g / l * sin(u[1])
  end

  G = 9.80665

  u0 = [theta0, dtheta0]  # initial angle and velocity
  tspan = (0.0, t)  # time span
  p = [G, l]  # parameters: gravity and length

  prob = ODEProblem(pendulum!, u0, tspan, p)
  sol = solve(prob, Tsit5(), reltol=1e-10)

  return sol[1, end]
end

@testset "nihalo1" begin
  l = 1.0
  t = 1.0
  theta0 = 0.1
  dtheta0 = 0.1
  n = 1000

  exact = exact_sol(l, t, theta0, dtheta0)
  @test isapprox(nihalo(l, t, theta0, dtheta0, n), exact, atol=1e-8)
end

@testset "nihalo_random" begin
  l = rand(1:10)
  t = rand(1:100)
  theta0 = rand(0:0.1:0.5)
  dtheta0 = rand(0:0.1:0.5)
  n = 10000

  exact = exact_sol(l, t, theta0, dtheta0)
  @test isapprox(nihalo(l, t, theta0, dtheta0, n), exact, atol=1e-4)
end
