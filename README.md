# Vaje pri predmetu Numerična matematika

## Navodila

To je projekt za delo pri predmetu [Numerična matematika](https://ucilnica.fri.uni-lj.si/course/view.php?id=117). Projekt je na začetku semestra prazen, tekom semestra bomo na vsakih vajah dodali rešitve. Gradiva z bolj podrobnim opisom so pripravljena v repozitoriju [z gradivi](https://nummat.gitlab.io/vaje-nummat/).

Predlagam, da si vsak naredi svoj Git repozitorij (na Gitlabu, Githubu ali kje drugje), kjer si uredi kodo in dokumentacijo, ki jo bomo pisali pri tem predmetu. Prav tako bomo git
uporabljali za oddajo domačih nalog.

## Navodila za hiter začetek

Ob začetku vaje si najprej ustvarimo direktorij oziroma paket za Julio[^1], kjer bo shranjeno naše delo

```shell
$ julia
julia> # pritisnemo ], da pridemo v način pkg
pkg> activate . # aktiviramo delovno okolje v trenutnem direktoriju (npr. vaje)
(vaje)pkg> generate VajaXY # generiramo ogrodje za novo vajo
(vaje)pkg> develop ./VajaXY  # dodamo pravkar generirano vajo v delovno okolje
(vaje) pkg> # pritisnemo tipko za brisanje nazaj, da zopet pridemo v navaden način
julia>
```

Zgornji ukazi ustvarijo direktorij `VajaXY` z osnovno struktura [paketa v Jiliji](https://pkgdocs.julialang.org/v1/creating-packages/). Za bolj obsežen projekt, lahko uporabite [šablono PkgTemplates](https://github.com/JuliaCI/PkgTemplates.jl) ali [PkgSkeleton.jl](https://github.com/tpapp/PkgSkeleton.jl).

```shell
julia> cd("VajaXY") # pritisnemo ;, da pridemo v način lupine
shell> tree .
.
├── Project.toml
└── src
    └── VajaXY.jl

1 directory, 2 files

```

Direktoriju dodamo še teste, skripte z demnostracijsko kodo in README dokument.

```shell
shell> mkdir test
shell> touch test/runtests.jl
shell> touch README.md
shell> mkdir scripts
shell> touch scripts/demo.jl
shell> tree .
.
├── Manifest.toml
├── Project.toml
├── scripts
│   └── demo.jl
├── src
│   └── VajaXY.jl
└── test
    └── runtests.jl
```

Ko je direktorij s kodo pripravljen lahko naložimo kodo iz `VajaXY.jl` v ukazni vrstici

```shell
julia> using VajaXY
julia> VajaXY.moja_super_funkcija()
```

Boljša možnost je, da kodo uporabimo v scripti npr. `scripts\demo.jl`.

```jl
# demo.jl vsebuje primere uporabe funkcije iz modula/paketa VajaXY
using VajaXY

VajaXY.moja_super_funkcija()
```

Scripto lahko poženemo z ukazom `ìnclude`.

```shell
julia> include("scripts/demo.jl")
```

Začetno strukturo paketa si lahko shranimo v šablono in uporabimo [PkgSkeleton.jl](https://github.com/tpapp/PkgSkeleton.jl) za generiranje novih paketov.

### Testi

Vstopna točka za teste je `test\runtests.jl`. Paket [Test](https://docs.julialang.org/en/v1/stdlib/Test/) omogoča pisanje enotskih testov, ki se lahko avtomatično izvedejo v sistemu [nenehne integracije (Continuous Integration)](https://en.wikipedia.org/wiki/Continuous_integration).

V juliji teste pišemo z makroji [@test](https://docs.julialang.org/en/v1/stdlib/Test/#Test.@test) in [@testset](https://docs.julialang.org/en/v1/stdlib/Test/#Test.@testset). Če `test/runtests.jl` lahko napišemo

```jl
using Test, VajaXY

@test VajaXY.funkcija_ki_vrne_ena() == 1
```

Lahko teste poženemo tako, da v `pkg` načinu poženemo ukaz `test VajeXY`

```shell
(vaje) pkg> test VajeXY

    Testing Running tests...
    Testing VajaXY tests passed
```

### Dokumentacija

Za pisanje dokumentacijo navadno uporabimo format [Markdown](https://en.wikipedia.org/wiki/Markdown). S paketom [Documenter](https://documenter.juliadocs.org/stable/) lahko komentarje v kodi in markdown dokumentente združimo in generiramo HTML ali PDF dokumentacijo s povezavo na izvorno kodo.

Za pripravo posameznih poročil lahko uporabite [IJulia](https://github.com/JuliaLang/IJulia.jl), [Weave.jl](https://github.com/JunoLab/Weave.jl), [Literate.jl](https://github.com/fredrikekre/Literate.jl) ali [Quadro](https://quarto.org/docs/computations/julia.html).

## Organizacija direktorijev

- `vaje` direktorij z vajami, ki služi kot [delovno okolje](https://pkgdocs.julialang.org/v1/environments/)
- `vaje/VajaXY` vsaka vaja ima svoj direktorij
- posamezen direktorij za vajo je organiziran kot paket s kodo, testi in dokumentacijo

        vaje
         └── Vaja01
           ├── Project.toml
           ├── README.md
           ├── src
           |   └─ Vaja01.jl
           ├── test
           |   └─ runtests.jl
           ├── doc
           |   ├─  makedocs.jl
           |   └─ index.md
           └─ scripts
              └─ demo.jl

## Delovno okolje

Za hitrejše in lažje delo z programskim jezikom `julia` uporabite [Revise](https://timholy.github.io/Revise.jl/stable/). Pred začetkom dela poženite

```julia
julia> using Revise
```
Zgornji ukaz bo poskrbel, da se bodo definicije funkcij sproti posodabljale, ko bomo spreminjali kodo v datotekah.

Če želite,  da se `Revise` zažene ob vsakem zagonu programa `julia`, lahko datoteko `startup.jl` namestite v `.julia/config/startup.jl`

```shell
$ cp startup.jl $HOME/.julia/config/startup.jl
```
### Generiranje PDF dokumentov

Za generiranje PDF dokumentov s paketi [Documenter](https://documenter.juliadocs.org/stable/) ali [Weave.jl](https://github.com/JunoLab/Weave.jl) je potrebno namestiti [TeX/LaTeX](https://tug.org/). Priporočam uporabo [TinyTeX](https://yihui.org/tinytex/).
Po [namestitvi](https://yihui.org/tinytex/#installation) tinytex, dodamo še nekaj `LaTeX` paketov, tako da v terminalu izvedemo naslednji ukaz

```
tlmgr install microtype upquote minted
```

## Povezave

- [Način dela za Gitlab (Gitlab Flow)](https://docs.gitlab.com/ee/topics/gitlab_flow.html).
- [Priporočila za stil Julia](https://docs.julialang.org/en/v1/manual/style-guide/).
- [Naveti za delo z Julijo](https://docs.julialang.org/en/v1/manual/workflow-tips/).

[^1]: Na vajah bomo uporabljali programski jezik Julia. Študentje lahko, če želijo, na vajah in za izdelavo domačih nalog uporabljajo katerikoli programski jezik.
