module DN2


export normal_CDF, bezier_curve_area, bezier

"""
  y = normal_CDF(x)

Izračunaj vrednost porazdelitvene funkcije Gaussove/normalne slučajne spremenljivke X ∼ N(0,1) za dano vrednost `x`.
"""
function normal_CDF(x, tol=10^-10)
  phi_prime(x) = exp(-0.5 * x^2) # e^(-x^2/2)
  MAXRANGE = 10

  if x > MAXRANGE
    return 1.0
  end
  if x < 0
    return 1 - normal_CDF(-x)
  end

  return adaptive_simpson(phi_prime, -MAXRANGE, x, tol) / koren(2 * pi())
end


"""
  y = koren(x)

Izračunaj vrednost kvadratnega korena danega števila `x˙. 
"""
function koren(a, tol=10^-10)
  x0 = zacetni_priblizek(a)
  for i = 1:100
    x = (x0 + a / x0) / 2
    if abs(x - x0) < tol * abs(x)
      return x
    end
    x0 = x
  end
  throw("Iteracija ne konvergira")
end


"""
  x0 = zacetni_priblizek(a)

Izračunaj začetni približek za kvadratni koren števila `a` z uporabo eksponenta števila s
plavajočo vejico. 
"""
function zacetni_priblizek(a)
  exp = div(exponent(a), 2) # sqrt(2^x) = 2^(x/2)

  if exp > 0
    return 1 << exp # magični začetni približek 2^(exp/2)
  else
    return 1 / (1 << -exp)
  end
end


"""
  p = pi()

Izračunaj približno vrednost konstante π.
"""
function pi()
  x = 2

  for i in 100:-1:1
    frac = i / (2 * i + 1)
    x = 2 + (frac * x)
  end

  return x
end


"""
  integral = adaptive_simpson(f, a, b, tolerance)

Izračunaj približno vrednost integrala funkcije `f` na intervalu `[a, b]` z dano toleranco `tolerance`
z uporabo adaptivne Simpsonove metode.
"""
function adaptive_simpson(f, a, b, tolerance)
  function simpson_rule(f, a, b)
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b))
  end

  function adaptive_simpson_recursive(f, a, b, tol, fa, fb, fMid, S)
    mid = (a + b) / 2
    d = (a + mid) / 2
    e = (mid + b) / 2
    fd = f(d)
    fe = f(e)
    Sleft = simpson_rule(f, a, mid)
    Sright = simpson_rule(f, mid, b)
    delta = abs(Sleft + Sright - S)

    if delta <= 15 * tol
      return Sleft + Sright + delta / 15
    end

    return adaptive_simpson_recursive(f, a, mid, tol / 2, fa, fMid, fd, Sleft) +
           adaptive_simpson_recursive(f, mid, b, tol / 2, fMid, fb, fe, Sright)
  end

  S = simpson_rule(f, a, b)
  return adaptive_simpson_recursive(f, a, b, tolerance, f(a), f(b), f((a + b) / 2), S)
end


"""
  area = bezier_curve_area()

Izračunaj ploščino zanke, ki jo vsebuje Bézierova krivulja dana s kontrolnim poligonom:
(0,0),(1,1),(2,3),(1,4),(0,4),(-1,3),(0,1),(1,0)
"""
function bezier_curve_area(tol=10^-10)
  control_points = [(0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (-1, 3), (0, 1), (1, 0)]
  t0, t1 = find_intersection_at_x(control_points, 0.5) # We know that the intersection is at x=0.5

  return 0.5 * adaptive_simpson(t -> bezier_integrand(control_points, t), t0, t1, tol)
end

"""
  b = binomial_coefficient(n, k)

Izračunaj binomski koeficient C(n, k) za dani `n` in `k`.
"""
function binomial_coefficient(n, k)
  if k < 0 || k > n
    return 0
  end
  if k == 0 || k == n
    return 1
  end

  return binomial_coefficient(n - 1, k - 1) + binomial_coefficient(n - 1, k)
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
"""
function bezier(control_points, t)
  n = length(control_points) - 1
  p = (0.0, 0.0)
  i = 0
  for control_point in control_points
    p = p .+ control_point .* bernstein_p(i, n, t)
    i += 1
  end

  return p
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
  t0, t1 = find_intersection_at_x(control_points, x_isect, tol)

Najdi parametra `t0` in `t1` pri katerih se Bézierova krivulja dana s kontrolnim poligonom `control_points`
seka v točki z x koordinato `x_isect`.
V primeru da presečišča ni, vrni `t0=0`, `t1=1`.
"""
function find_intersection_at_x(control_points, x_isect, tol=0.001)
  t0 = 0.0
  t1 = 1.0
  y_isect = nothing

  for t in range(0, stop=1, length=1000)
    x, y = bezier(control_points, t)
    if abs(x - x_isect) < tol
      if isnothing(y_isect)
        y_isect = y
        t0 = t
      elseif abs(y_isect - y) < tol
        # NOTE: Only the first point with x=x_isect is considered
        t1 = t
        break
      end
    end
  end

  return t0, t1
end


end # module DN2