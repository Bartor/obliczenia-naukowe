# Autor: Bartosz Rajczyk

# Każda funkcja w tym programie bierze za argument jakiś typ danych liczb zmiennoprzecinkowych
# Tablica "types" trzyma wykorzysytwane przez nas typy danych
types = [Float16, Float32, Float64]

function machEps(type)
    number = one(type)
    while one(number) + number / 2 != one(number)
        number /= 2
    end
    number
end

println("Machine Epsilon")
foreach(type -> println("f: ", machEps(type), " eps: ", eps(type)), types)

function eta(type)
    number = one(type)
    while number / 2 != zero(type)
        number /= 2
    end
    number
end

println("Eta")
foreach(type -> println("f: ", eta(type), " nextfloat: ", nextfloat(zero(type))), types)

function maxNumber(type)
    number = one(type)
    while !isinf(number * 2)
        number *= 2
    end
    # aktualnie number = inf/2, musimy uzupełnić brakującą wartość
    gap = number / 2
    while !isinf(number + gap) && gap >= one(type)
        number += gap
        gap /= 2
    end
    number
end

println("Max")
foreach(type -> println("f: ", maxNumber(type), " floatmax: ", floatmax(type)), types)