module Vaje02

import Base: getindex, setindex!, size
import Base: *, \

#' # Definirajmo podatkovno strukturo za tridiagonalne matrike
"""
Podatkovna strukutra za tridiagonalno matriko.
"""
struct Tridiag
  sp  # Vektor elementov na spodnji diagonali
  d   # Vektor elementov na diagonali
  zg  # Vektor elementov na zgornji diagonalie
end

export Tridiag

"""
Vrni dimenzije tridiagonalne matrike.
"""
size(T::Tridiag) = (length(T.d), length(T.d))

"""
  elt = T[i, j]

Vrni element v `i`-ti vrstici in `j`-tem stoplcu tridiagonalne matrike `T`.
"""
function getindex(T::Tridiag, i, j)
  n, _m = size(T)
  if (i < 1) || (i > n) || (j < 1) || (j > n)
    throw("Indeksi so izven obsega matrike")
  end

  if i == j - 1
    return T.zg[i]
  elseif i == j
    return T.d[i]
  elseif i == j + 1
    return T.sp[j]
  else
    zero(T.d[1])
  end
end

"""
Nastavi vrednost elementa v `i`-ti vrstici in `j`-tem stoplcu tridiagonalne matrike `T` na `x`.

"""
function setindex!(T::Tridiag, x, i, j)
  n, _m = size(T)
  if (i < 1) || (i > n) || (j < 1) || (j > n)
    throw("Indeksi so izven obsega matrike")
  end

  if i == j - 1
    T.zg[i] = x
  elseif i == j
    T.d[i] = x
  elseif i == j + 1
    T.sp[j] = x
  else
    throw("Indeksi so izven tridiagonalne matrike")
  end
end

"""
Izracunaj produkt tridiagonalne matrike z vektorjem.
"""
function *(T::Tridiag, x::Vector)
  y = zero(x)
  n = length(T.d)

  y[1] = T[1, 1] * x[1] + T[1, 2] * x[2]
  for i = 2:n-1
    y[i] = T[i, i-1] * x[i-1] + T[i, i] * x[i] + T[i, i+1] * x[i+1]
  end
  y[n] = T[n, n-1] * x[n-1] + T[n, n] * x[n]

  return y
end


function \(T::Tridiag, b::Vector)
  # i-ti korak algoritma
  # l = a[i][i-1] / a[i-1][i-1]
  # a[i][i] -> a[i][i] - l * a[i][i-1]
  # b[i] -> b[i] * l * b[i-1]

  n, _m = size(T)
  T = Tridiag(T.sp, deepcopy(T.d), T.zg)
  b = float(deepcopy(b))
  # eliminacija
  for i = 2:n
    l = T[i, i-1] / T[i-1, i-1]
    T[i, i] = T[i, i] - l * T[i-1, i]
    b[i] = b[i] - l * b[i-1]
  end

  # obratna substitucija
  x = float(zero(b))
  x[n] = b[n] / T[n, n]
  for i = (n-1):-1:1
    x[i] = (b[i] - T[i, i+1] * x[i+1]) / T[i, i]
  end

  return x
end

export getindex
export *
export \


end # module Vaje02
