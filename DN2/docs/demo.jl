#' # DN2 - Porazdelitvena funkcija normalne slučajne spremenljivke, Ploščina Bézierove krivulje; Matija Ojo; 63200205

#' ## Porazdelitvena funkcija normalne slučajne spremenljivke 

#' Porazdelitvena funkcija normalne slučajne spremenljivke je definirana kot:

#' $$F(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} e^{-\frac{t^2}{2}} dt$$.

#' V praksi se uporablja za izračun verjetnosti, da bo vrednost slučajne spremenljivke padla v določen interval.

#' Za izračun porazdelitvene funkcije normalne slučajne spremenljivke se uporablja le numerične metode, saj integrala ne moremo izračunati analitično.
#' V nalogi je uporabljena [sestavljenova Simpsonova 1/3 metoda](https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_1/3_rule) za izračun integrala.
#' Sestavljena Simpsonova metoda je osnovana na [Simpsonovem pravilu](https://en.wikipedia.org/wiki/Simpson%27s_rule), ki numerično izračuna vrednost integrala tako da evalvira vrednost funckije na vmesnih točkah.
#' Velja: $\int_a^b f(x) \approx \frac{b - a}{6} (f(a) + 4f(\frac{a+b}{2}) + f(b))$.
#' Sestavljeno Simpsonovo 1/3 pravilo pa interval $[a, b]$ razdeli na $n$ podintervalov in izračuna integral na vsakem podintervalu.
#' Velja namreč: $\int_a^b f(x) \approx \frac{1}{3} h (f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + 2f(x_4) + \dots + 2f(x_{n-2}))$.

#' ### Rezultati

#' Poglejmo si graf relativne napake med porazdelitveno funkcijo normalne slučajne spremenljivke izračunano
#' s sestavljeno Simpsonovo metodo ter z metodo izračuna integrala z uporabo Gaussove kvadrature iz paketa
#' `FastGaussQuadrature` in jo smatramo za "točno" vrednost.

using DN2
using Plots
using FastGaussQuadrature

function accurate_gauss_CDF(x)
  f(x) = exp(-x^2 / 2) / sqrt(2 * pi)
  a = -10
  b = x
  n = 60
  xs, ws = gausslegendre(n)
  xs = (b - a) / 2 * xs .+ (b + a) / 2
  ws = (b - a) / 2 * ws
  return sum(ws .* f.(xs))
end

xs = range(-1, stop=1, length=1000)
ys = [gaussian_CDF(x) for x in xs]
true_ys = [accurate_gauss_CDF(x) for x in xs]
diffs = [abs(ys[i] - true_ys[i]) / abs(true_ys[i]) for i in 1:length(xs)]
plot(xs, diffs, label="Relativna napaka")



#' ## Ploščina Bézierove krivulje

#' Bézierovo krivuljo omejujeo kontrolne točke $P_i$ in je definirana kot:

#' $$B(t) = \sum_{i=0}^{n} b_{i,n}(t) P_i; \quad 0 \leq t \leq 1$$

#' kjer je:

#' $$b_{i,n}(t) = {{n}\choose{i}} t^i (1 - t)^{n-i}; \quad i = 0 \ldots 1$$

#' i-ti Bernsteinov bazni polinom stopnje $n$.

#' V nalogi je Bézierova krivulja podana z naslendjimi kontrolnimi točkami:

control_points = [(0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (-1, 3), (0, 1), (1, 0)];

#' Poglejmo si njen graf

t = range(0, stop=1, length=1000)
x_values = [bezier(control_points, t_value)[1] for t_value in t]
y_values = [bezier(control_points, t_value)[2] for t_value in t]
plot(x_values, y_values, xlabel="x", ylabel="y", title="Bézierova krivulja z zanko", legend=:none)


#' Ploščino zanke se izračuna tako da se najde samo-presek krivulje (začetek in konec zanke) in meji
#' uporabimo pri sestavljeni Simpsonovi metodi za izračun ploščine.

#' Meji so izračunani s pomočjo Newtonove metode z relativno natančnostjo $10^{-10}$.

#' ## Rezultati 

#' Ploščina omenjene zanke je enaka:

bezier_curve_area()

#'