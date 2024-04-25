module Vaje07

using LinearAlgebra

export newton


"""
  x, it = newton(fdf, x0)

Poisci resitev nelinearnega sistema enacb z Newtonovo metodo.

## Argumenti

- fdf ... funkcija, ki racuna vrednosti funkcije in Jacobijeve matrike za dani argument
- x0 ... zacetni priblizek
- maxit ... maksimalno stevilo iteracij
- tol ... natancnost, s katero racunamo priblizek (razlika dveh sosednjih priblizkov)
- eps ... natancnost za vrednost f(x)
"""
function newton(fdf, x0, maxit=100, tol=1e-10, eps=1e-10)
  for i = 1:maxit
    f, Jf = fdf(x0)
    x = x0 - Jf \ f

    if norm(f, Inf) < eps || norm(x - x0, Inf) < tol
      return x, i
    end

    x0 = x
  end
  throw("Newtonova metoda ne konvergira po $maxit korakih.")
end

end # module Vaje07
