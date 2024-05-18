using Vaje09
using Plots


#' # Numericna integracija

#' Splosna oblika numericne interacije je:

#' $\int_a^b f(x) dx \approx \sum_{i=1}^n w_i f(x_i)$

#' kjer so $x_i$ vozlišča in $w_i$ uteži.

#' Integracijsko pravilo zahteva da dolocimo utezi ter vozlisca.

#' Uporabili bomo metode:

#' - trapezno pravilo
#' - Simpsonovo pravilo
#' - Gauss-Legendrove kvadrature

#' Ocenili bomo napako in primerjali rezultate.

#' ## Sestavljeno trapezno pravilo

#' Osnovno trapzeno pravilo:

#' $\int_a^b f(x) dx \approx \frac{b - a}{2} (f(a) + f(b)) + R_f$

#' kjer je $R_f$ napaka, odvisna od $f$.

#' Sestavljeno trapzeno pravilo z $n$ podintervali:

#' $\int_a^b f(x) dx \approx \sum_{i=1}^n \frac{h}{2} (f(x_i) + f(x_{i+1})) + R_f = \frac{h}{2} (f(a) + 2f(a + h) + \dots + 2f(b - h) + f(b)) + R_f$

#' ### Izračunajmo integral sin na [0, 1]

#' $\int_0^1 \sin(x) dx = -\cos(1) + \cos(0) = 1 - \cos(1)$

Ip = 1 - cos(1) # prava vrednost

#'

#' Izračunajmo integral s trapezno metodo z $n = 5$ koraki.
It5 = integral(sin, trapez, (0, 1, 5))

#'

relative_err5 = (It5 - Ip) / Ip

#'

#' Izračunajmo integral s trapezno metodo z $n = 10$ koraki.
It10 = integral(sin, trapez, (0, 1, 10))

#'

relative_err10 = (It10 - Ip) / Ip

#'


#' #### Graf napake v odvisnosti od števila korakov trapezne metode.

napaka_t(n) = Ip - integral(sin, trapez, (0, 1, n))

tabn = 2 .^ (1:10)

scatter(tabn .+ 1, abs.(napaka_t.(tabn)), xaxis=:log10, yaxis=:log10, label="trapez")

#'

# Približen red metode
red_trapez = hcat(ones(length(tabn)), log10.(tabn)) \ log10.(abs.(napaka_t.(tabn)))

#'

#' ## Simpsonovo pravilo

#' #### Graf napake za Simpsonovo metodo

napaka_s(n) = Ip - integral(sin, simpson, (0, 1, n))

scatter!(2 * tabn .+ 1, abs.(napaka_s.(tabn)), xaxis=:log10, yaxis=:log10, label="Simpson")

#'

#' Približen red Simpsonove metode
red_simpson = hcat(ones(length(tabn)), log10.(2 * tabn .+ 1)) \ log10.(abs.(napaka_s.(tabn)))

#'

#' ## Gauss-Legendre
Igl5 = integral(sin, gl, (0, 1, 5))
gl_err5 = Igl5 - Ip

#'

#' Graf napake za Gauss-Legendrove kvadrature

napaka_gl(n) = Ip - integral(sin, gl, (0, 1, n))

scatter!(tabn, abs.(napaka_gl.(tabn)) .+ 5e-16, xaxis=:log10, yaxis=:log10,
  label="Gauss-Legendre")

#'
