#' # Potencna metoda za iskanje lastnih vrednosti

#' A matrika dimenzij n x n
#' Iscemo lastne pare lambda, v:
#' $A\vec{v} = λ \vec{v}, \vec{v} \neq 0$

#' Potencna metoda:
#' $x^{(n+1)} = \frac{A x^{(n)}}{|| A x^{(n)} ||}$

#' Zaporedje $x^{(n)}$ konvergira k lastnemu vektorju za lastno vrednost, ki je po absolutni vrednosti najvecja.

#' Invariantna mera Markovske verige

#' Markovska veriga:
#' $X_1, X_2, \ldots, X_n, \ldots$
#' $P(X_n=i | X_{n-1}=j) = p_{ji}$

#' Matrika prehodnih verjetnosti:

#' $P = [p_{11}, p_{12}, \ldots, p_{1n}; p_{21} \dots; p_{n1} \dots p_{nn}]$

#' $\sum_{j=1}^{n} p_{ij} = 1$

#' $P * \vec{1} = \vec{1} \Rightarrow P \vec{1} = λ \vec{1}, λ=1$

#' $X_1 ∼ (p_1^0 p_1^0 \dots p_n^0)$

#' $\vec{p}_0$ ... zacetna porazdelitev za $X_1$

#' $\vec{p}_1 = P^T\vec{p}_0$ ... porazdelitev za $X_2$

#' $\vec{p}_i^1 = P(X_2 = i) = \sum_{j=1} P(X_2 = i | X_1 = j) * P(X_1 = j) = \sum p_{ji} * p_j^0$

#' $lim_{n \to \infty} (P^T)^n * \vec{p}^0 = \vec{p}$

#' Kako bi izracunali lambdo ce poznamo (priblizek) za lastni vekor?

#' $x, A * x \approx λ * x$

#' $x = [x_1, x_2, \ldots, x_n]$

#' $Ax = [λ x_1, λ x_2, \ldots, λ x_n]$

#' Imamo $n$ enacb s katerimi lahko izracunamo lastno vrednost.
#' Izberemo maksimalni element (po absolutni vrednosti) vektor $Ax$ delimo s tem elementom.

#' $\lambda = \frac{(Ax)_i}{x_i}$ za $i$ pri katerem je maksimalni element.

#' Izboljsava: vektor $x$ delimo z najvecjim elementom po abs. vrednosti.
#' Tako bo maksimalni element v $x$ enak $1$, oziroma v $Ax$ bo maksimalni element enak $\lambda$.

#' max(\frac{x}{max(x_i)})=1$

using Vaje05
using LinearAlgebra
using Plots

P = [0.1 0.4 0.5; 0 0.5 0.5; 0.2 0 0.8]

#' Invariantna mera je lastni vektor za $P^T$ za lastno vrednost $1$.

v, lambda = potencna(P', ones(3), 100, 1e-10)

#' Uporabljeno v Google PageRank algoritmu, kjer so definirali Markovsko verigo, ki predstavlja gibanje po povezavah med spletnimi stranmi.

#' Invariatna mera za Markovsko matriko na 6 vozliscih.


bipartite_M = [
    0 0.3 0 0.4 0 0.3;
    0.1 0 0.2 0 0.7 0;
    0 0.5 0 0.2 0 0.3;
    0.4 0 0.2 0 0.4 0;
    0 0.5 0 0.2 0 0.3;
    0.4 0 0.2 0 0.4 0
]


eigen(bipartite_M')

#'

v, lambda = potencna(bipartite_M', rand(6), 100, 1e-10)

# scatter(v, title="Invariantna porazdelitev za MV na 6 tockah")

#' Ce ima matrika A vec razlicnih vrednosti, ki so po absolutni vrednosti najvecje,
#' potem potencna metoda ne konvergira za vsak zacetni priblizek

#' Konkretno: potencna metoda je imela tezave pri konvergenci, pri lastni vrednostih $1$ in $-1$.

#' Resitev: premik: $A - \delta I$

#' Ce so $\lambda_1, \dots, \lambda_n$ lastne vrednosti matrike A, potem so
#' $\lambda_1 - \delta, \dots, \lambda_n - \delta$ lastne vrednosti matrike $A - \delta I$.

B = bipartite_M' + I

#'

v, lambda = potencna(B', rand(6))

#' 

scatter(v, title="Invariantna porazdelitev za MV na 6 tockah")

#'