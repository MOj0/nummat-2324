module Vaje03

using SparseArrays

export laplaceova_matrika, desne_strani, RobniProblemPravokotnikLaplace
export resi, resi_iter

"""
    L = laplaceova_matrika(n, m)

Zapiši matriko za Laplaceov operator v 2D na pravokotnem območju. 
Matrika `L` je matrika sistema enačb za diskretizirano laplaceovo enačbo

```math
u_{i-1,j}+u_{i,j-1} - 4u_{ij} + u_{i+1,j}+u_{i,j+1} = 0.
```
"""
function laplaceova_matrika(n, m)
  L = spdiagm(0 => -4 * ones(n), 1 => ones(n - 1), -1 => ones(n - 1))
  I = spdiagm(0 => ones(n)) # identiteta
  A = spzeros(n * m, n * m)
  for j = 1:m
    k = ((j-1)*n+1):(j*n) # indeksi v j-tem bloku  
    A[k, k] = L
    if j < m
      A[k, k.+n] = I
      A[k.+n, k] = I
    end
  end
  return A
end


"""
  b = desne_strani(s, z, l, d)

Generiraj desne strani za robni problem na pravokotniku [a, b]x[c, d].
s ... vektor robnih pogojev spodaj (y = c)
z ... vektor robnih pogojev zgoraj (y = d)
l ... vektor r. p. levo (x = a)
d ... vektor r. p. desno (x = b)
"""
function desne_strani(s, z, l, d)
  n = length(s)
  m = length(l)
  b = zeros(n * m)
  b[1:n] -= s # j = 1
  b[n:n:end] -= d # i = n
  b[end-n+1:end] -= z # j = m
  b[1:n:end-n+1] -= l # i = 1
  return b
end

struct RobniProblemPravokotnikLaplace
  meje # meje [a, b] x [c, d] v vektorju [a, b, c, d]
  robni_pogoji # funkcije na robu [fs, fz, fl, fd] f(a, y) = fl(y), f(x, c) = fs(x), ...
end

"""
  Z = resi(rp, h)

Poišči približno rešitev RP za Laplaceovo enačbo za pravokotnik, 
tako da pravokotnik razdelimo na mrežo s korakom h.
"""
function resi(rp::RobniProblemPravokotnikLaplace, h)
  a, b, c, d = rp.meje
  fs, fz, fl, fd = rp.robni_pogoji
  n = Integer(floor((b - a) / h) - 1)
  m = Integer(floor((d - c) / h) - 1)
  # Laplaceova matrika
  L = laplaceova_matrika(m, n)
  # desne strani
  x = range(a, b, n + 2)
  y = range(c, d, m + 2)
  Z = zeros(m + 2, n + 2)
  Z[end, :] = fs.(x)
  Z[1, :] = fz.(x)
  Z[:, 1] = fl.(y)
  Z[:, end] = fd.(y)
  b = desne_strani(Z[end, 2:end-1], Z[1, 2:end-1], Z[2:end-1, 1], Z[2:end-1, end])
  # resi sistem
  x = L \ b
  # resitev preoblikuj v matriko
  Z[2:end-1, 2:end-1] = reshape(x, m, n)
  return Z
end

"""
Reši robni problem z iterativno metodo SOR
"""
function resi_iter(rp::RobniProblemPravokotnikLaplace, h, w)
  a, b, c, d = rp.meje
  fs, fz, fl, fd = rp.robni_pogoji
  n = Integer(floor((b - a) / h) - 1)
  m = Integer(floor((d - c) / h) - 1)
  # robni pogoji
  x = range(a, b, n + 2)
  y = range(c, d, m + 2)
  Z = zeros(m + 2, n + 2)
  Z[end, :] = fs.(x)
  Z[1, :] = fz.(x)
  Z[:, 1] = fl.(y)
  Z[:, end] = fd.(y)
  iteracija(Z, w)
end
"""
Poišči rešitev Laplacoeve enačbe s SOR iteracijo.
"""
function iteracija(Z, w, maxit=1000)
  n, m = size(Z)
  Z0 = copy(Z)
  for it = 1:maxit
    for i = 2:n-1
      for j = 2:m-1
        # Gauss-Seidel
        Z[i, j] = (Z[i+1, j] + Z[i, j+1] + Z[i-1, j] + Z[i, j-1]) / 4
        Z[i, j] = w * Z[i, j] + (1 - w) * Z0[i, j] # SOR popravek
      end
    end
    if isapprox(Z, Z0, atol=1e-5)
      println("Število iteracij: $it")
      return Z
    end
    copyto!(Z0, Z)
  end
  throw("Iteracija ni konvergirala v $maxit korakih")
end
end # module Vaje03
