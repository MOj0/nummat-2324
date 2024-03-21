#' # Vaja 03 - Minimalne ploskve

using Vaje03
using LinearAlgebra
using Plots


#' Poišči minimalno ploskev napeto na rob s pravokotnim tlorisom.

#' ## Laplaceova matrika

#' 1. Ploskev predstavi z vrednostmi `z` koordinate na kvadratni mreži.
#' 2. Matriko, ki predstavlja ploskev spremeni v vektor. 
#' 3. Zapiši linearen sistem za Laplaceovo enačbo za minimalno ploskev.
#' 4. Reši sistem. Uporabi razpršene matrike (paket `SparseArrays`).
#' 5. Grafično predstavi rešitev (ploskev)
#' 6. Grafično predstavi, kako se ob LU razcepu razpršena Laplaceova matrika napolni. 


#' Funkcijo $f(x, y)$ obravnavamo v vozliscih mreze

#' $k(i, j) = j + m * (i - 1)$
#' $1, 1 -> 1$


#' Diskretizacija Laplaceove enacbe (koncne diference)

#' $f'(x) \approx  \frac{f(x, + h/2) - f(x - h/2)}{h}$

#' $f''(x) \approx \frac{\frac{f(x+h) - f(x)}{h} - \frac{f(x) - f(x-h)}{h}}{h} = \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}$


#' Robni problem (diferencialni operator)

#' $△ f = 0$

#' $[a, b] \times [c, d] \dots \text{obmocje}$

#' Robni pogoji
#' $f(a, y) = f_1(y)$
#' $f(b, y) = f_2(y)$
#' $f(x, c) = f_3(x)$
#' $f(x, d) = f_4(x)$


#' Podatki za robni problem
#' zapakiramo v podatkovno strukturo: RobniProblemPravokotnik
#' Funkcija, ki resi robni problem: resi(robni_problem, podatki)

rp = RobniProblemPravokotnikLaplace([0, pi, -1, 1], [sin, sin, x -> 0.0, x -> 0.0])
Z = resi(rp, 0.1)
surface(Z[2:end-1, 2:end-1])

#' Izris ploskve


#' ### Napolnitev Laplaceove matrike

L = laplaceova_matrika(10, 10)
spy(L)

#' Oblika L


F = lu(L)
spy(F.L)

#' Faktorja L in U v razcepu imata več neničelnih elementov, kot originalna matrika.


spy(F.U)
#' Zgornje trikotna matrika


#' # Iterativne metode

Z = resi_iter(rp, 0.1, 1) # Gauss-Seidlova iteracija
surface(Z)

#' Slika ploskve pri $w=1$

Z = resi_iter(rp, 0.1, 1.85)
surface(Z)