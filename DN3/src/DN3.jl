module DN3

export nihalo, nihalo_izracun, nihalo_harmonicno

"""
  y = rk4(f, y0, t0, h)

Izračunaj naslednjo vrednost `y` z Runge-Kutta metodo 4. reda.
"""
function rk4(f, y0, h)
  k1 = h * f(y0)
  k2 = h * f(y0 + k1 / 2)
  k3 = h * f(y0 + k2 / 2)
  k4 = h * f(y0 + k3)
  return y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6
end

""""
  y = nihalo(l, t, theta0, dtheta0, n)

Izračunaj končno stanje nihala ob času `t`, z uporabo Runge-Kutta metode 4. reda.
Nihalo je dolžine `l`, začetni odmik (odmik ob času 0) je `theta0` in začetna kotna hitrost (𝜃′(0)) je `dtheta0`.
Interval [0, t] je razdeljen na `n` podintervalov enake dolžine.
"""
function nihalo_izracun(l, t, theta0, dtheta0, n)
  G = 9.80665
  h = t / n
  y = [theta0, dtheta0]

  # theta' = y
  # y' = -g/l * sin(theta)
  f(x) = [x[2], -G / l * sin(x[1])]

  curr_t = 0
  for _ in 1:n
    y = rk4(f, y, h)
    curr_t += h
  end

  return y
end


""""
  odmik = nihalo(l, t, theta0, dtheta0, n)

Izračunaj odmik nihala ob času `t`.
Parametri:
- `l`: dolžina nihala
- `t`: čas
- `theta0`: začetni odmik
- `dtheta0`: začetna kotna hitrost
- `n`: število korakov
"""
function nihalo(l, t, theta0, dtheta0, n)
  return nihalo_izracun(l, t, theta0, dtheta0, n)[1] # Vrne samo odmik
end


"""
  odmik = nihalo_harmonicno(l, t, theta0, n)

Izračunaj odmik nihala ob času `t`, z uporabo harmoničnega (sinusnega) nihanja.
Nihalo je dolžine `l`, začetni odmik (odmik ob času 0) je `theta0`.
Interval [0, t] je razdeljen na `n` podintervalov enake dolžine.
"""
function nihalo_harmonicno(l, t, theta0, n)
  G = 9.80665
  h = t / n
  curr_t = 0
  odmik = theta0

  for _ in 1:n
    odmik = theta0 * cos(curr_t * sqrt(G / l))
    curr_t += h
  end

  return odmik
end


end # module DN3
