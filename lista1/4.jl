# Autor: Bartosz Rajczyk

"""
Funkcja wyznaczająca taką liczbę zmiennoprzecinkową 1 < x < 2, że 1 * 1/x != 1
"""
function a()
    current = one(Float64)
    while nextfloat(current) * (one(Float64)/nextfloat(current)) == one(Float64) && current < 2
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = a()
println(bitstring(result))
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))

"""
Funkcja wyznaczająca najmniejszą taką liczbę zmiennoprzecinkową 0 < x, że 1 * 1/x != 1
"""
function b()
    current = nextfloat(zero(Float64))
    while nextfloat(current) * (one(Float64)/nextfloat(current)) == one(Float64)
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = b()
println(bitstring(result))
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))
