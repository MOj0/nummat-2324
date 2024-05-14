module DN2


export gaussian_CDF, bezier_curve_area, bezier

"""
  y = gaussian_CDF(x)

Izračunaj vrednost porazdelitvene funkcije Gaussove/normalne slučajne spremenljivke X ∼ N(0,1) za dano vrednost `x`.
"""
function gaussian_CDF(x)
  phi_prime(x) = exp(-0.5 * x^2)
  return 1 / sqrt(2 * pi) * simpson_rule(phi_prime, -10, x, 600)
end


"""
  s = simpson_rule(f, a, b, n)

Izračunaj vrednost integrala funkcije `f` na intervalu `[a, b]` s sestavljenim
Simpsonovim 1/3 pravilom na `n` točkah (`n` je sodo).
"""
function simpson_rule(f, a, b, n)
  h = (b - a) / (2 * n)
  s = f(a) + f(b) + 4 * f(a + h)
  for i in 1:n-1
    s += 2 * f(a + 2 * i * h) + 4 * f(a + 2 * i * h + h)
  end
  s *= h / 3
  return s
end


"""
  area = bezier_curve_area()

Izračunaj ploščino zanke, ki jo vsebuje Bézierova krivulja dana s kontrolnim poligonom:
(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)
"""
function bezier_curve_area(tol=10^-10)
  control_points = [(0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (-1, 3), (0, 1), (1, 0)]
  t0, t1 = find_intersection_at_x(control_points, 0.5, 1.0)

  return 0.5 * simpson_rule(t -> bezier_integrand(control_points, t), t0, t1, 500)
end

"""
  b = binomial_coefficient(n, k)

Izračunaj binomski koeficient C(n, k) za dani `n` in `k`.
"""
function binomial_coefficient(n, k)
  if (k > n - k)
    k = n - k
  end
  res = 1
  for i in 0:k-1
    res *= (n - i) // (i + 1)
  end
  return res
end


"""
  p = bernstein_p(i, n, t)

Izračunaj vrednost i-tega Bernsteinovega baznega polinoma stopnje `n` pri vrednosti `t`.
"""
function bernstein_p(i, n, t)
  binomial_coefficient(n, i) * t^i * (1 - t)^(n - i)
end

"""
  p = bezier(control_points, t)

Izračunaj vrednost Bézierove krivulje dane s kontrolnim poligonom `control_points` pri vrednosti `t`.
Uporabi de Casteljaujev algoritem.
"""
function bezier(control_points, t)
  n = length(control_points)
  points = [(float(control_points[i][1]), float(control_points[i][2])) for i in 1:n]

  for i in 1:n-1
    for j in 1:n-i
      points[j] = (1 - t) .* points[j] .+ t .* points[j+1]
    end
  end

  return points[1]
end

"""
  d = bezier_derivative(control_points, t)

Izračunaj odvod Bézierove krivulje dane s kontrolnim poligonom `control_points` pri vrednosti `t`.
"""
function bezier_derivative(control_points, t)
  n = length(control_points) - 1
  d = (0.0, 0.0)
  i = 0
  for (p, pNext) in zip(control_points[1:end-1], control_points[2:end])
    d = d .+ (pNext .- p) .* bernstein_p(i, n - 1, t)
    i += 1
  end

  return d .* n
end

"""
  b = bezier_integrand(control_points, t)

Izračunaj integrand Bézierove krivulje dane s kontrolnim poligonom `control_points` pri vrednosti `t`.
Uporablja se formula za ploščino krivočrtnega trikotnika pod krivuljo, ki predpostavi da sta končni točki enaki:
bezier(t_0) == bezier(t_end).
"""
function bezier_integrand(control_points, t)
  x, y = bezier(control_points, t)
  x_prime, y_prime = bezier_derivative(control_points, t)
  return x * y_prime - y * x_prime
end


"""
  t0, t1 = find_intersection_at_x(control_points, x_isect, y_isect, tol)

Najdi parametra `t0` in `t1` pri katerih se Bézierova krivulja dana s kontrolnim poligonom `control_points`
seka v točki z x koordinato `x_isect` in y koordinato manjso od `y_upper`.
"""
function find_intersection_at_x(control_points, x_isect, y_upper, tol=10^-10)
  t0 = 0.0
  t1 = 1.0

  for t in range(0, stop=1, length=10)
    x, y = bezier(control_points, t)
    if abs(x - x_isect) < 0.3 && y < y_upper
      if t0 == 0.0  # 0.0 can be represented exactly with floats, so this is fine
        t0 = t
      else
        t1 = t
        break
      end
    end
  end

  maxiter = 100
  for _ in 1:maxiter
    x0, y0 = bezier(control_points, t0)
    dx0, _ = bezier_derivative(control_points, t0)
    if abs(x0 - x_isect) / abs(x_isect) < tol && y0 < y_upper
      break
    end
    t0 -= (x0 - x_isect) / dx0
  end

  for _ in 1:maxiter
    x1, y1 = bezier(control_points, t1)
    dx1, _ = bezier_derivative(control_points, t1)
    if abs(x1 - x_isect) / abs(x_isect) < tol && y1 < y_upper
      break
    end
    t1 -= (x1 - x_isect) / dx1
  end

  return t0, t1
end


end # module DN2