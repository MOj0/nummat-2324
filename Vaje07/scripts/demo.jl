#' # Nelinearni sistemi enacb

#' $n \times n$ sistem

#' $f_1(x_1, x_2, \dots, x_n) = 0$

#' $f_2(x_1, x_2, \dots, x_n) = 0$

#' $\vdots$

#' $f_n(x_1, x_2, \dots, x_n) = 0$

#' $\vec{F}(\vec{x}) = \vec{0}$

#' ## Newtonova metoda

#' Za $n=1$

#' $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$

#' Za $n>1$

#' Kaj je odvod $\vec{F}(\vec{x})$?

#' $\vec{F}(\vec{x}): \mathbb{R}^n \rightarrow \mathbb{R}^n$?


#' Jacobijeva matrika

#' $J\vec{F}: \mathbb{R}^n \rightarrow \mathbb{R}^{n \times n}$

#' $$J\vec{F} = D\vec{F} = \partial \vec{F} = \begin{bmatrix}
#' \frac{\partial u_1}{\partial x_1} & 
#' \frac{\partial u_1}{\partial x_2} & 
#' \frac{\partial u_1}{\partial x_3}
#' \frac{\partial u_2}{\partial x_1} & 
#' \frac{\partial u_2}{\partial x_2} & 
#' \frac{\partial u_2}{\partial x_3} \\[1ex]
#' \frac{\partial u_3}{\partial x_1} & 
#' \frac{\partial u_3}{\partial x_2} & 
#' \frac{\partial u_3}{\partial x_3}
#'  \end{bmatrix}$$

#' Torej je enacba: $\vec{x}^{(n+1} = \vec{x}^{n} - J \vec{F}(\vec{x})^{-1} * \vec{F}(\vec{x}^{(n)})$

#' ## Primer: Presecisce 2 kroznic

#' $(x+1)^2 + y^2 = 4$

#' $(x-1)^2 + y^2 = 4$

#' Preuredimo enacbi v vektorski zapis:

#' $$\vec{F}(\vec{x}) = \begin{bmatrix}
#' (x_1+1)^2 + x_2^2 - 4 \\
#' (x_1-1)^2 + x_2^2 - 4
#' \end{bmatrix}$$


#' Izracunamo Jacobijevo matriko:


#' $$J\vec{F} = \begin{bmatrix}
#' 2(x_1+1) & 2x_2 \\
#' 2(x_1-1) & 2x_2
#' \end{bmatrix}$$


#' Kako vemo, da smo nasli resitev?

#' $\vec{x}'$ priblizek, $\vec{x}$ pa prava resitev.

#' 1) $\vec{F}(\vec{x}') = 0$
#' $ ||F(\vec{x}')|| < \epsilon$

#' 2) $||\vec{x}' - \vec{x}|| < \delta$
#' $||\vec{x}^{(n+1)} - \vec{x}^{(n)}|| < \delta$

#' Ustrezno moramo kombinirati oba pogoja ("po x" ali "po y")


#' ## Primeri

#' - presecisce dveh parametricnih krivulj
#' - presecisce krivulje s ploskvijo
#' - minimum funkcije vec spremenljivk (razdalja med dvema krivuljama / krivuljo in ploskvijo)


#' Naredimo primer za razdaljo med dvema krivuljama

#' Za dve parametricno podani krivulji poisci minimalno razdaljo med dvema tockama na krivuljah.

#' Razdalja: $d_2(t, s) = (x-1(t) - x_2(s))^2 + (y_1(t) - y_2(s))^2$

#' Lokalni ekstremi so v stacionarnih tockah funkcije $d2$.

#' $\frac{d_2(t, s)}{\partial t} = 0$

#' $\frac{d_2(t, s)}{\partial s} = 0$

#' Sistem je torej: $grad(d_2)(\vec{ts}) = \vec{0}$

#' $J(grad(d_2)) = H(d_2)$

#' Vzamimo 2 elipsi


using Vaje07
using Plots
using ForwardDiff
using LinearAlgebra

tocka(a) = tuple(a...)

K1(t) = [2 * cos(t) + 1 / 3, sin(t) + 0.25]
K2(s) = [cos(s) / 3 - sin(s) / 2, cos(s) / 3 + sin(s) / 2]
t = LinRange(0, 2 * pi, 60);
plot(tocka.(K1.(t)), label="K1")
plot!(tocka.(K2.(t)), label="K2")

#' Iscemo mimimum kvadrata razdalje

function d2(ts)
  delta = K1(ts[1]) - K2(ts[2])
  return dot(delta, delta)
end

# ForwardDiff.gradient(d2, [1, 2])

f(ts) = ForwardDiff.gradient(d2, ts)
JF(ts) = ForwardDiff.hessian(d2, ts)
fdf(ts) = (f(ts), JF(ts))


ts, x = newton(fdf, [0, 0])

#'

scatter!(tocka(K1(ts[1])), label="K1")
scatter!(tocka(K2(ts[1])), label="K2")

#' Graf funkcije razdalje

contour(range(-pi, 2pi, 100), range(-pi, 2pi, 100), (t, s) -> d2([t, s]))

scatter!(tuple(ts...))

#' Drug zacenti prilblizek...

ts, x = newton(fdf, [2, 1])
scatter!(tuple(ts...))