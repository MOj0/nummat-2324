module DN2


export normal_CDF

"""
  y = normal_CDF(x)

Izračunaj vrednost Gaussove/normalne porazdelitvene funkcije X ∼ N(0,1) za dano vrednost `x`.
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


end # module DN2


