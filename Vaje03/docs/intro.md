

n * m blocna matrika, ki ima bloke velikosti m*m.


Funkcijo f(x, y) obravnavamo v vozliscih mreze

f(x_i, y_j) = z_{ij} = ?


k(i, j) = j + m * (i - 1)

1, 1 -> 1


Diskretizacija Laplaceove enacbe (koncne diference)

f'(x) \approx  \frac{f(x, + h/2) - f(x - h/2)}{h}

f''(x) \approx \frac{\frac{f(x+h) - f(x)}{h} - \frac{f(x) - f(x-h)}{h}}{h} = \frac{f(x+h) - 2f(x) + f(x-h)}{h^2}