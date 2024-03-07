using Vaje02
using LinearAlgebra
using Plots


#' # Vaja2 - linearni sistemi

#' Resi sistem linearnih enacb
#' ```
#' x + 2y - 7 = 1
#' 2x - y - 3z = 2
#' x - z = 3
#' ```
#' Sistem najprej pretvorimo v matricno obliko
#' $Ax = b$.

A = [1 2 -1; 2 -1 -3; 1 0 -1]  # matrika sistema
b = [1, 2, 3] # desne strani

#' resitev sistema $A x = b; (A \textbackslash b = A^{-1} b)$
x = A \ b

#' Preizkus naredimo tako, da preverimo, ce je $Ax - b$ enak 0.

zero = A * x - b


#' Sistem lahko resimo tudi z LU razcepom, tako da sistem $L U x = b$ prevedemo na 2 sistema
#' $L y = b$ in $Ux = y$

#' L - spodnje trikotna matrika

#' U - zgornje trikotna matrika

#' $A x = L U x=b;Ux = y; Ly = b; x= U \textbackslash (L \textbackslash b)$ 


L, U, p = lu(A)

x = U \ (L \ b[p])

zero2 = norm(A * x - b, Inf)

#' Julia faktorje razcepa zapakira v en objekt.

F = lu(A)
x = F \ b

#' ## Tridiagonalne matrike

T = Tridiag([1, 2], [3, 4, 5], [6, 7])

elt = T[1, 1] + T[2, 2]
println(T[1, 1], " + ", T[2, 2], " = ", elt)

#' Produkt matrike z vektorjem

T * [1, 2, 3]


L = Tridiag([1, 1], [-2, -2, -2], [1, 1])
b1 = [1, 1, 1, 1]
x = L \ b1
# println(x)



#' ## Slucajni sprehod

#' $X_n = \sum_{i=1}^{n} Bin(p)$

#' Naloga: Izracunaj pricakovano stevilo korakov ki ga potrebuje slucajni prehod da doseze n ali -n.

"Generator nakljucnih stevil porazdeljenih kot Bin(p)"
rand_bin(p) = (rand() < p) ? 1 : -1

[rand_bin(0.6) for i = 1:10]

sprehod(p, n) = cumsum([rand_bin(p) for i = 1:n])

# scatter(sprehod(0.5, 1000), label="Sprehod 1")

#' Markovska veriga: prehod v naslednjo stanje je odvisno samo od prejsnjega

#' Markovske verige se da predstaviti z matriko prehodov.

#' Markovska veriga, predstavljena z matriko za slucajni sprehod:

# [1 0   ...   0]
# [q 0 p  ...  0]
# [  q 0 p ... 0]
#      ...
# [            p]
# [          q 0]

# Sistem, ki ga moramo resit

# [1-p       ] [k_1     [1
# [-q 1 -p   ]  .        .
# [  ...     ]  .    =   .
# [        -p]  .        .
# [      -q 1] k_n]      1]


n = 100
p = 0.5
T = Tridiag(-p * ones(2 * n - 2), ones(2 * n - 1), -(1 - p) * ones(2 * n - 2))
k = T \ ones(2 * n - 1)
scatter(-n:n, k)