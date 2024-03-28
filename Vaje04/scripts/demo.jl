#' # Aproksimacija z linearnim modelom

using Plots
using LinearAlgebra

data = readlines(open("num_data.txt", "r"))
data = strip.(data)

filter!(l -> l[1] != '#', data)
data = [split(line, r"\s+") for line in data]
data = [[parse(Float64, x) for x in line] for line in data]

filter!(l -> l[5] > 0, data)
t = [l[4] for l in data]
co2 = [l[5] for l in data]
scatter(t, co2, title="Atmosferski CO2 na Mauna Loa",
    xlabel="leto", ylabel="parts per milion (ppm)", label="co2",
    markersize=1)

#' 

t0 = 2000

A = hcat(ones(size(t)), t .- t0, (t .- t0) .^ 2, cos.(2pi * t), sin.(2pi * t))

p = A \ co2

plot(t, p[2] .+ 2 * p[3] * (t .- t0), title="Letne spremembe CO2")

#' 

model(t) = p[1] + p[2] * (t - t0) + p[3] * (t - t0)^2 + p[4] * cos(2 * pi * t) + p[5] * sin(2 * pi * t)
print(model.([2020, 2030, 2050]))
