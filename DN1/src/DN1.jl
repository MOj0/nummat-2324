module DN1

using Plots

export Zlepek, interpoliraj, vrednost, izirsi_zlepek

"""
Podatkovna struktura za zlepek.
Vsebuje koeifiente za vsak odsek in točke, kjer je zlepek interpoliran.
"""
struct Zlepek
  coefficients
  interpolation_points_x
end


"""
  Z = interpoliraj(x, y)

  Izračunaj koeficiente zlepka na interpolacijskih točkah in vrni strukturo `Zlepek`
"""
function interpoliraj(x, y)
  n = length(x)
  a = copy(y)

  b = zeros(n - 1)
  d = zeros(n - 1)
  h = [x[i+1] - x[i] for i in 1:n-1]

  alpha = zeros(n - 1)
  for i in 2:n-1
    alpha[i] = 3 / h[i] * (a[i+1] - a[i]) - 3 / h[i-1] * (a[i] - a[i-1])
  end

  c = zeros(n)
  l = zeros(n)
  mu = zeros(n)
  z = zeros(n)

  l[1] = 1.0

  for i in 2:n-1
    l[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]
    mu[i] = h[i] / l[i]
    z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]
  end

  l[n] = 1
  z[n] = 0
  c[n] = 0

  for j in n-1:-1:1
    c[j] = z[j] - mu[j] * c[j+1]
    b[j] = (a[j+1] - a[j]) / h[j] - (h[j] * (c[j+1] + 2 * c[j])) / 3
    d[j] = (c[j+1] - c[j]) / (3 * h[j])
  end

  coefficients = zeros(n - 1, 4)
  for i in 1:n-1
    coefficients[i, 1] = a[i]
    coefficients[i, 2] = b[i]
    coefficients[i, 3] = c[i]
    coefficients[i, 4] = d[i]
  end

  return Zlepek(coefficients, x)
end

"""
  y = vrednost(spline, x)

  Izracunaj vrednost zlepka pri tocki x.
"""
function vrednost(spline::Zlepek, x::Real)
  i = length(spline.interpolation_points_x) - 1
  for index in 1:length(spline.interpolation_points_x)
    if x < spline.interpolation_points_x[index]
      i = index - 1
      break
    end
  end

  a = spline.coefficients[i, 1]
  b = spline.coefficients[i, 2]
  c = spline.coefficients[i, 3]
  d = spline.coefficients[i, 4]
  x0 = spline.interpolation_points_x[i]

  return a + b * (x - x0) + c * (x - x0)^2 + d * (x - x0)^3
end

"""
  plot_spline(Z)

  Izriši zlepek.
"""
function izirsi_zlepek(Z::Zlepek)
  n = length(Z.interpolation_points_x)
  colors = ["blue", "red"]

  plot(title="Zlepek na $n točkah", label="")

  for i in 1:length(Z.interpolation_points_x)-1
    xs = range(Z.interpolation_points_x[i], stop=Z.interpolation_points_x[i+1], length=100)
    ys = [vrednost(Z, x) for x in xs]

    plot!(xs, ys, label="", linewidth=2, color=colors[i%2+1])
  end

  x_orig = Z.interpolation_points_x
  y_orig = [vrednost(Z, x) for x in x_orig]
  plot!(x_orig, y_orig, seriestype=:scatter, label="")
end


end # module DN1
