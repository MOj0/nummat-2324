module Vaje05

"""
  v, lambda = potencna(A, x0, maxit, tol)

Izracunaj lastni vektor `v` in lastno vrednost `lambda` matrike `A` z uporabo potencne metode
z zacetnim priblizkom za lastni vektor `x0`.
"""
function potencna(A, x0, maxit=100, tol=1e-10)
  x = x0
  _, max_idx = findmax(abs, x)
  ls = x[max_idx]
  x = x / ls

  for i = 1:maxit
    x = A * x
    _, maxj = findmax(abs, x)
    ln = x[maxj]
    x = x / ln  # normiramo, maksimalni element v x je 1

    if abs(ln - ls) < tol
      println("Konvergirano po $i iteracijah.")
      return (x, ln)
    end

    ls = ln
  end

  throw("Potencna metoda ni konvergirala v $maxit korakih.")
end

export potencna

end # module Vaje05
