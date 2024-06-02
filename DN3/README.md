# Avtor: Matija Ojo

## Matematično nihalo

### Opis naloge:

Kotni odmik 𝜃(𝑡) (v radianih) pri nedušenem nihanju nitnega nihala opišemo z diferencialno enačbo

```math
g/l * sin(𝜃(𝑡)) + 𝜃''(𝑡) = 0, 𝜃(0)=𝜃_0, 𝜃′(0)=𝜃′_0,
``` 

kjer je 𝑔=9.80665𝑚/𝑠^2 težni pospešek in 𝑙 dolžina nihala.
Napišite funkcijo nihalo, ki računa odmik nihala ob določenem času.
Enačbo drugega reda prevedite na sistem prvega reda in računajte z metodo Runge-Kutta četrtega reda:

```math
k1 = h f(x_n, y_n)
k2 = h f(x_n + h/2, y_n + k1 / 2)
k3 = h f(x_n + h/2, y_n + k2 / 2)
k4 = h f(x_n + h, y_n + k3)
y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4) / 6
``` 


Klic funkcije naj bo oblike odmik=nihalo(l,t,theta0,dtheta0,n) kjer je

- `odmik` enak odmiku nihala ob času `t`,
- dolžina nihala je `l`,
- začetni odmik (odmik ob času 0)) je `theta0`
- in začetna kotna hitrost (𝜃′(0)) je `dtheta0`,
- interval [0,𝑡] razdelimo na `n` podintervalov enake dolžine.

Primerjajte rešitev z nihanjem harmoničnega nihala. Za razliko od harmoničnega nihala (sinusno nihanje), je pri matematičnem nihalu nihajni čas odvisen od začetnih pogojev (energije). Narišite graf, ki predstavlja, kako se nihajni čas spreminja z energijo nihala.


## Uporaba kode

Implementirana sta funkciji `nihalo(l,t,theta0,dtheta0,n)`, ki z uporabo metode Runge-Kutta 4. reda izračuna odmik matematičnega nihala v času `t`,
ter `nihalo_harmonicno(l, t, theta0, n)`, ki izračuna odmik harmoničnega nihala v času `t`.


## Poganjanje testov

```shell
nummat-2324$ julia
julia> # pritisnemo ], da pridemo v način pkg
pkg> activate DN3
(DN1)pkg> test
```

## Ustvarjanje poročila

```shell
nummat-2324$ julia --project=. DN3/docs/make.jl
```