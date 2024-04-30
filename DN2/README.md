# Avtor: Matija Ojo

## Porazdelitvena funkcija normalne slučajne spremenljivke

### Opis naloge:

Napišite učinkovito funkcijo, ki izračuna vrednosti porazdelitvene funkcije za standardno normalno porazdeljeno slučajno spremenljivko 𝑋∼𝑁(0,1).

## Uporaba kode

Implementirana je funkcija `normal_CDF(x, tol=10^-10)`, ki preko adaptivnega Simpsonovega pravila izračuna vrednost porazdelitvene funkcije normalne slučajne spremenljivke.


## Ploščina Bézierove krivulje

### Opis naloge:

Izračunajte ploščino zanke, ki jo omejuje Bézierova krivulja dana s kontrolnim poligonom:

(0,0),(1,1),(2,3),(1,4),(0,4),(−1,3),(0,1),(1,0).

## Uporaba kode

Implementirana je funkcija `normal_CDF(tol=10^-10)`, ki najde samo-presečišče v krivulji (začetek in konec zanke) in preko adaptivnega Simpsonovega pravila izračuna površino znotraj zanke.


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