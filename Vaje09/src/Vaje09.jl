module Vaje09

using FastGaussQuadrature

"""
Poišči Gauss-Legendrove kvadrature z `n` točkami za interval [a, b]
"""
function gl(a, b, n)
  x, w = gausslegendre(n)
  x = (b - a) / 2 * x .+ (a + b) / 2
  w = (b - a) / 2 * w
  return x, w
end

"""
    x, w = trapez(a, b, n)

Izračunaj vozlišča `x` in uteži `w` za sestavljeno trapenzno pravilo na 
intervalu `[a, b]` z `n` podintervali.
"""
function trapez(a, b, n)
  x = range(a, b, n + 1)
  h = x[2] - x[1]
  w = ones(n + 1) * h
  w[1] /= 2
  w[end] /= 2
  return x, w
end

"""
    x, w = simpson(a, b, n)

Izračunaj vozlišča `x` in uteži `w` za sestavljeno Simpsonovo pravilo na
intervalu `[a, b]` z `n` podintervali.
"""
function simpson(a, b, n)
  x = range(a, b, 2n + 1)
  h = x[2] - x[1]
  w = h / 3 * ones(2n + 1)
  w[2:2:end-1] *= 4
  w[3:2:end-1] *= 2
  return x, w
end

"""
    I = integral(fun, metoda, parametri)

Izračunaj integral funkcije `fun` z metodo `metoda` s `parametri`.
Argument `metoda` je funkcija, ki vrne vektorja vozlišč in uteži.
"""
function integral(fun, metoda, parametri)
  x, w = metoda(parametri...)
  return w' * fun.(x)
end

export trapez, simpson, gl, integral
end # module Vaje09
