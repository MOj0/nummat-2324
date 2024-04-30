#' # DN2 - Porazdelitvena funkcija normalne slučajne spremenljivke, Ploščina Bézierove krivulje; Matija Ojo; 63200205

#' ## Porazdelitvena funkcija normalne slučajne spremenljivke 

#' Porazdelitvena funkcija normalne slučajne spremenljivke je definirana kot:

#' $$F(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} e^{-\frac{t^2}{2}} dt$$.

#' V praksi se uporablja za izračun verjetnosti, da bo vrednost slučajne spremenljivke padla v določen interval.

#' Za izračun porazdelitvene funkcije normalne slučajne spremenljivke se uporablja le numerične metode, saj integrala ne moremo izračunati analitično.
#' V nalogi je uporabljena [adaptivna Simpsonova metoda](https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method) za izračun integrala.
#' Adaptivna Simpsonova metoda je osnovana na [Simpsonovem pravilu](https://en.wikipedia.org/wiki/Simpson%27s_rule), ki numerično izračuna vrednost integrala tako da evalvira vrednost funckije na vmesnih točkah.
#' Velja: $\int_a^b = \frac{b - a}{6} (f(a) + 4f(\frac{a+b}{2}) + f(b))$.
#' Adaptivna metoda pa razdeli interval na dva dela in izračuna integral na vsakem delu posebej, in njuno vsoto uporabi kot naslendji približek.
#' Če je razlika med dobljenim pribiližkom integrala in prejšnjim približkom večja od določene tolerance, se interval razdeli na še manjše dele, sicer vrne izračunani približek.

#' ### Rezultati

#' Poglejmo si graf porazdelitvene funkcije normalne slučajne spremenljivke za različne vrednosti $x$,
#' pri čemer je vsaka vrednost funckije izračunana z adaptivno Simpsonovo metodo.

using DN2
using Plots

xs = range(-4, stop=4, length=1000)
ys = [normal_CDF(x) for x in xs]
plot(xs, ys, xlabel="x", ylabel="F(x)", title="Porazdelitvena funkcija N(0, 1)", legend=:none)


#' ## Ploščina Bézierove krivulje

#' Bézierovo krivuljo omejujeo kontrolne točke $P_i$ in je definirana kot:

#' $$B(t) = \sum_{i=0}^{n} b_{i,n}(t) P_i; \quad 0 \leq t \leq 1$$

#' kjer je:

#' $$b_{i,n}(t) = {{n}\choose{i}} t^i (1 - t)^{n-i}; \quad i = 0 \ldots 1$$

#' i-ti Bernsteinov bazni polinom stopnje $n$.

#' V nalogi je Bézierova krivulja podana z naslendjimi kontrolnimi točkami:

control_points = [(0, 0), (1, 1), (2, 3), (1, 4), (0, 4), (-1, 3), (0, 1), (1, 0)]

#' Poglejmo si njen graf

t = range(0, stop=1, length=1000)
x_values = [bezier(control_points, t_value)[1] for t_value in t]
y_values = [bezier(control_points, t_value)[2] for t_value in t]
plot(x_values, y_values, xlabel="x", ylabel="y", title="Bézierova krivulja z zanko", legend=:none)


#' Da izračunamo ploščino zanke, moramo najprej najti presečišče, kjer se zanka začne (in tudi konča).
#' Opazimo lahko da ima presečišče x-koordinato enako $0.5$.
#' To lahko uporabimo pri  iskanju vrednosti neodvisne spremenljivke `t` tega presečišča.
#' S preprostim algoritmom lahko ugotovimo da se pri `t0=0.075` zanka začne in se konča pri `t1=0.924`.
#' Ti dve meji uporabimo pri izračunu ploščine zanke, natančneje, pri adaptivni Simpsonovi metodi.


#' ## Rezultati 

#' Ploščina omenjene zanke je enaka:

bezier_curve_area()

#'