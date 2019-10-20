# Autor: Bartosz Rajczyk

# Każda funkcja w tym programie bierze za argument jakiś typ danych liczb zmiennoprzecinkowych
# Tablica "types" trzyma wykorzysytwane przez nas typy danych
types = [Float16, Float32, Float64]

function funnyMachEps(type)
    typeOne = one(type)
    3*typeOne*((4*typeOne)/(3*typeOne) - typeOne) - typeOne
end

foreach(type -> println("f: ", funnyMachEps(type), " eps: ", eps(type)), types)
