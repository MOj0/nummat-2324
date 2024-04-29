# TODO: Spremeni vse!


# Naravni zlepek

## Avtor: Matija Ojo

## Opis naloge:

Imamo $n$ interpolacijskih točk $(x_i, f_i),\quad i=1,2\ \ldots\  n$.
Naravni interpolacijski kubični zlepek $S$ je funkcija, ki izpolnjuje naslednje pogoje:
1. graf zlepka gre skozi interpolacijske točke: $S(x_i) = f_i,\quad i=1,2\ \ldots n$,
2. je polinom stopnje 3 ali manj na vsakem podintervalu $[x_i, x_{i+1}],\quad i=1,2\ \ldots n-1$,
3. je dvakrat zvezno odvedljiva funkcija na interpolacijskem intervalu $[x_1,x_n]$,
4. drugi odvod v začetni in končni točki je enak 0: $S''(x_1)=S''(x_n)=0$.


Napišite funkcijo `Z = interpoliraj(x, y)`, ki izračuna koeficient polinoma $S_i$ in vrne rezultat tipa `Zlepek`.

Napišite funkcijo `y = vrednost(Z, x)`, ki vrne vrednost zlepka v dani točki $x$.

Napišite funkcijo `izirsi_zlepek(Z)`, ki nariše graf zlepka, tako da različne odseke izmenično nariše z rdečo in modro barvo.

## Uporaba kode

Zlepek je možno izračunati tako, da se definira seznam stevil `x` in `y`, ki morata biti enake dolzine,
nato pa se klice funkcijo `interpoliraj(x, y)`.

Ta funkcija vrne objekt tipa `Zlepek`, s katerim je mozno klicati funkcijo `y = vrednost(Z, x)`, ki izracuna
vrednost zlepka v tocki `x`.

Zlepek je mozno tudi izrisati s klicem funkcije `izirsi_zlepek(Z)`.


## Poganjanje testov

```shell
nummat-2324$ julia
julia> # pritisnemo ], da pridemo v način pkg
pkg> activate DN2
(DN1)pkg> test
```

## Ustvarjanje poročila

```shell
nummat-2324$ julia --project=. DN2/docs/make.jl
```