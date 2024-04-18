using Vaje06
using Random, Plots, LinearAlgebra

#' # Inverzna iteracija

A = [3 1 1; 1 1 3; 1 3 4]
F = lu(A)
resi(b) = F \ b
v, lambda = inviter(resi, 3)

#' Preizkusimo
norm(A * v - lambda * v, Inf)

#'

eigen(A)


#' # Spektralno grucenje

#' Imamo oblak tock v $R^n$, ki jih preslikamo v drug prostor, tako da so podobne tocke blizje.

#' Podobnostni graf: "podobne" tocke so povezane.

#' Graf $\epsilon$ okolic ($G_\epsilon$): tocki sta povezani, ce sta oddaljeni manj kot $\epsilon > 0$.

#' Za vsak graf definiramo Laplaceovo matriko, za matriko sosednosti A:

#' $L = D - A =  \begin{pmatrix} st(1) & 0 & \dots & -1 & \dots & -1 \\ 0 & st(2) & \dots & -1   \end{pmatrix}$

#' Matrika $D$ je diagonalna, kjer je $d_ii = st(i)$ stopnja vozlisca $i$.

#' Za nove koordinate tock uporabimo komponente lastnih vektorjev za L.

#' Graf podamo na vec nacinov:

#' - seznam povezav
#' - seznam sosedov
#' - matrika sosednosti: $a_ij = 1 \iff (i, j) \in E, a_ij = 0 \text{sicer}$


#' Kot primer vzamimo mnozico treh normalno porazdeljenih oblakov tock v ravnimi.


n = 100

Random.seed!(12)

g1 = randn(2, n) .+ (1, -2)
g2 = randn(2, n) .+ (-3, -1)
g3 = randn(2, n) .+ (0, 1)

scatter(g1[1, :], g1[2, :], label="1. gruca")
scatter!(g2[1, :], g2[2, :], label="2. gruca")
scatter!(g3[1, :], g2[2, :], label="3. gruca")

#' Najprej zgradimo podobnostni matriki

oblak = hcat(g1, g2, g3)
eps = 1
A = graf_eps(oblak, eps);

#' in izracunamo Laplaceovo matriko grafa:
L = laplace(A);

#' Poiscimo lastne vektorje matrike L za majhne lastne vrednosti.

F = lu(L + 1e-5 * I)
v, lambda = inviterqr(b -> F \ b, 300, 10, 1000)

lambda

#' Za koordinate uporabimo lastne vektorje Laplaceove matrike

a = 2
b = a + 1
scatter(v[1:100, a], v[1:100, b], label="1. gruca")
scatter!(v[101:200, a], v[101:200, b], label="2. gruca")
scatter!(v[201:300, a], v[201:300, b], label="3. gruca")

#' Vgrajena funkcija eigen

v, l = eigen(Matrix(L))