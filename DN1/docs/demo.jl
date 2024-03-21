#' # DN1 - Naravni zlepek

#' Imamo $n$ interpolacijskih točk $(x_i, f_i)$, $i=1,2,...,n$.
#' Naravni interpolacijski kubični zlepek $S$ je funkcija, ki izpolnjuje naslednje pogoje:
#' 1. $S(x_i) = f_i, i=1,2,\dots,n$
#' 2. je polinom stopnje 3 ali manj na vsakem podintervalu $[x_i, x_{i+1}], i=1,2,\dots,n-1$
#' 3. je dvakrat zvezno odvedljiva funkcija na interpolacijskem intervalu $[x_1,x_n]$
#' 4. $S''(x_1)=S''(x_n)=0$


#' ## Izračun koeficientov

#' Za izračun koeficientov obstaja več načinov.
#' Koeficiente $a_i, b_i, c_i, d_i$ kubičnega zlepka $S$ lahko izračunamo tako, da uporabimo tridiagonalno matriko $A$ in vektor $d$.
#' Vendar pa bomo v nalogi uporabili algoritem, ki izračuna omenjene koeficiente.
#' Opis algorima samega presega obseg te naloge, zato ga ne bomo podrobneje opisovali.
#' Zavedati se moramo le, da algoritem zagotovi naslednje lastnosti zlepka (poleg že omenjenih na začetku):
#' 1. $S_i(x_i) = f_i = S_{i-1}(x_i)$
#' 2. $S'_i(x_i) = S'_{i-1}(x_i)$
#' 3. $S''_i(x_i) = S''_{i-1}(x_i)$

#' Radoveden bralec pa lahko o algoritmu več izve na [Wikipediji](https://en.wikipedia.org/wiki/Spline_(mathematics)#Algorithm_for_computing_natural_cubic_splines).

#' ## Rezultati

using DN1

#' Definirajmo nekaj poljubnih točk
x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
y = [1.0, 3.0, 1.0, 2.0, 0.0, 6.0]
s = interpoliraj(x, y)
p1 = vrednost(s, 5)
println("p1=", p1)
#' Veljati mora $p1=6.0$, saj je točka $(5.0, 6.0)$ bila podana kot interpolacijska točka.

p2 = vrednost(s, 4.99)
println("p2=", p2)
#' Veljati mora $p2≈6.0$, saj je točka z $x=4.99$ blizu točki $(5.0, 6.0)$.

izirsi_zlepek(s)

#' Kot vidimo iz slike, se zlepek prilega interpolacijskim točkam.


#' Poglejmo si še en primer.

x = range(0, stop=2 * pi, length=5)
y = sin.(x)
s = interpoliraj(x, y)

izirsi_zlepek(s)

#' Kot vidimo iz slike, je zlepek uspešno interpoliral funkcijo sinus.