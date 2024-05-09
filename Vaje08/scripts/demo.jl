#' # Interpolacija

#' $(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)$


#' f je interpolant na $(x_i, y_i)$, ce $f(x_i) = y_i$ za vse $i$.

#' ## Cebiševe tocke

#' Za njih velja da je interpolacija numericno stabilna.

#' ## Zlepki

#' Funkcije z vec predpisi


#' # Hermitov kubicni zlepek

#' Poleg vrednosti funkcije, je na voljo tudi vrednost odvoda.

#' $(x_1, y_1, dy_1), (x_2, y_2, dy_2), \ldots, (x_n, y_n, dy_n)$

#' So zlepki C1, saj drugi odvodi niso zvezni.

#' Podatkov: $2n$, $n-1$ intervalov.

#' Na vsakem intervalu imamo 4 podatke.

#' Enacbe so med intervali med seboj $\textbf{neodvisne}$.

#' Ko resujemo sistem enacb s $4$ neznankami: $y_0, y_1, dy_0, dy_1$ racunamo naslednji interpolacijski polinom:

#' $p_3(x) = a + bx + cx^2 + dx^3$

#' Najprej resimo problem na $[0, 1]$.

#' Imamo naslednjo Hermitovo bazo polinomov.

#' $h_{00}(t)(p(0)) = 1$
#' $h_{00}(t)(p(1)) = 0$
#' $h_{00}(t)(p'(0)) = 0$
#' $h_{00}(t)(p'(1)) = 0$

#' $h_{01}(t)(p(0)) = 0$
#' $h_{01}(t)(p(1)) = 1$
#' $h_{01}(t)(p'(0)) = 0$
#' $h_{01}(t)(p'(1)) = 0$

#' $h_{10}(t)(p(0)) = 0$
#' $h_{10}(t)(p(1)) = 0$
#' $h_{10}(t)(p'(0)) = 1$
#' $h_{10}(t)(p'(1)) = 0$

#' $h_{11}(t)(p(0)) = 0$
#' $h_{11}(t)(p(1)) = 0$
#' $h_{11}(t)(p'(0)) = 0$
#' $h_{11}(t)(p'(1)) = 1$

#' $h_{00}:$

#' $h_{ij}(0) = a_{ij}$
#' $h_{ij}(1) = a_{ij} + b_{ij} + c_{ij} + d_{ij}$
#' $h'_{ij}(0) = b_{ij}$
#' $h'_{ij}(1) = b_{ij} + 2 c_{ij} + 3d_{ij}$



#' $$A = \begin{bmatrix}
#' 1 & 0 & 0 & 0 \\
#' 1 & 1 & 1 & 1 \\
#' 0 & 1 & 0 & 0 \\
#' 0 & 1 & 2 &  3 
#'  \end{bmatrix} * 
#' \begin{bmatrix}
#' 1 & 0 & 0 & 0 \\
#' 0 & 1 & 0 & 0 \\
#' 0 & 0 & 1 & 0 \\
#' 0 & 0 & 0 &  1 
#'  \end{bmatrix}$$



#' $$A^{-1} = \begin{bmatrix}
#' 1 & 0 & 0 & 0 \\
#' 0 & 0 & 1 & 0 \\
#' -3 & 3 & -2 & -1 \\
#' 2 & -2 & 1 &  1 
#'  \end{bmatrix}$$


#' Preslikati moramo interval $[0, 1] \leftrightarrow [x_i, x_{i+1}]$

#' $t: 0$
#' $f(x(t)) = y_i$
#' $\frac{d}{dt} = dy_i$

#' $t: 1$
#' $f(x(t)) = y_{i+1}$
#' $\frac{d}{dt} = (x_{i+1} - x_i)$


#' $x: x_i$
#' $f(x) = y_i$
#' $f'(x) = dy_i$

#' $x: x_{i+1}$
#' $f(x) = y_{i+1}$
#' $f'(x) = dy_{i+1}$

#' Ce t definiramo kot:
#' $t = \frac{x - x_i}{x_{i+1} - x_i}$

#' Potem so vrednosti v robnih tockah:

#' $t(x_i) = 0$
#' $t(x_{i+1}) = 1$

#' Odvodi so naslednji:
#' v $t$ odvajamo po $x$ in dobimo:

#' $\frac{dt}{dx} = \frac{1}{x_{i+1} - x_i}$

#' Inverz je reciprocna vrednost:

#' $\frac{dx}{dt} = \frac{x_{i+1} - x_i}{1}$

#' $\frac{d}{dt} f(x(t)) = f'(x) * x'(t) = f'(x) * (x_{i+1} - x_i)$



using Vaje08
using Plots

x = range(0, 5pi, 7)
y = sin.(x)
dy = cos.(x)

scatter(x, y, label="Podatki")
Z = HermitovZlepek(x, y, dy)
plot!(x -> Z(x), 0, 5pi, label="Hermitov zlepek")
plot!(sin, 0, 5pi, label="sin")