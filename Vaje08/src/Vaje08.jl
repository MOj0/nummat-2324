module Vaje08

export hermiteInt, HermitovZlepek, vrednost

# bazne funkcije na [0, 1]
h00(t) = 1 + t^2 * (-3 + 2 * t)
h01(t) = t^2 * (3 - 2 * t)
h10(t) = t * (1 + t * (-2 + t))
h11(t) = t^2 * (-1 + t)


"""
  y = hermiteInt(x, xint, y, dy)

Izracunaj vrednost Hermitovega kubicnega polinoma p v tocki `x`, ki interpolira
podatke `p(xint[1] == y[1])`, `p(xint[2] == y[2])` in
`p'(xint[1]) == dy[1]`, `p'(xint[2]) == dy[2]`.
"""
function hermiteInt(x, xint, y, dy)
  dx = xint[2] - xint[1]
  t = (x - xint[1]) / dx
  return y[1] * h00(t) + y[2] * h01(t) + dx * (dy[1] * h10(t) + dy[2] * h11(t))
end

"""
Podatkovna struktura, ki hrani podatke za Hermitov kubicni zlepek v interpolacijskih tockah `x`
z danimi vrednostmi `y` in vrednostjo odvoda `dy`.
"""
struct HermitovZlepek
  x
  y
  dy
end

(Z::HermitovZlepek)(x) = vrednost(x, Z)

"""
  z = vrednost(x, Z)

Izracunaj vrednost Hermitovega kubicnega zlepka `Z` v dani toki `x`.
"""
function vrednost(x, Z::HermitovZlepek)
  if x == first(Z.x)
    return first(Z.y)
  end
  if x == last(Z.x)
    return last(Z.y)
  end

  i = searchsortedfirst(Z.x, x)
  if i > lastindex(Z.x) || i == firstindex(Z.x)
    throw(BoundsError(Z, i))
  end

  return hermiteInt(x, Z.x[i-1:i], Z.y[i-1:i], Z.dy[i-1:i])
end


end # module Vaje08
