# Autor: Bartosz Rajczyk

"""
Funkcja wyznaczająca taką liczbę zmiennoprzecinkową 1 < x < 2, że 1 * 1/x != 1

type - typ danych, dla którego ją wyznaczyć
"""
function a(type)
    current = one(type)
    while nextfloat(current) * (one(type)/nextfloat(current)) == one(type)
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = a(Float64)
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))

"""
Funkcja wyznaczająca najmniejszą taką liczbę zmiennoprzecinkową 0 < x, że 1 * 1/x != 1

type - typ danych, dla którego ją wyznaczyć
"""
function b(type)
    current = nextfloat(zero(type))
    while nextfloat(current) * (one(type)/nextfloat(current)) == one(type)
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = b(Float64)
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))
