# Autor: Bartosz Rajczyk

"""
Funkcja testująca manualnie, czy każda liczba pomiędzy podanymi dwoma
znajduje się w podanej odległości między sobą.

first - liczba, od której zacząć sprawdzanie
last - liczba, na której zakończyć sprawdzanie
step - tetsowana odległość pomiędzy kolejnymi liczbami
"""
function experimentalTest(first :: Float64, last :: Float64, step :: Float64)
    current = nextfloat(first)
    while (current < last)
        if (nextfloat(current) != current + step)
            return false;
        end
        current = nextfloat(current)
    end
    true
end

# poprawne, ale zdecydowanie za długie
# println(experimentalTest(1.0, 2, 2^(-52)))


"""
Funkcja testująca na podstawie znajomości reprezentacji liczb zmiennoprzecinkowych
w IEEE754 podwójnej precyzji, czy wszystkie kolejne liczby w podanym przedziale są w podanej
odległości od siebie.

first - liczba, od której zacząć sprawdzanie
last - liczba, na której zakończyć sprawdzanie
step - tetsowana odległość pomiędzy kolejnymi liczbami

"""
function bitstringTest(first :: Float64, last :: Float64, step :: Float64)
    last = prevfloat(last)
    firstExp = SubString(bitstring(first), 2:12)
    lastExp = SubString(bitstring(last), 2:12)
    # jeżeli eksponenty nie są równe, nie będziemy mieli równego rozmieszczenia
    if (firstExp != lastExp)
        return false
    end

    # eksponenta to zapis z biasem 1023 w przypadku float64, a zapis mantysy ma 52 miejsca, więc do 2^(-52) dokładności
    if (((2.0^(parse(Int, firstExp, base = 2) - 1023))*2.0^(-52)) != step)
        return false
    end

    return true
end

println(bitstringTest(1.0, 2.0, 2^(-52)))

"""
Funckaj wyznaczająca na podstawie znajomości zapisu liczb zmiennoprzecinkowych
w IEEE754 podwójnej precyzji odległość między kolejnymi liczbami w zakresie.

first - początek zakresu
last - koniec zakresu
"""
function numberSpread(first :: Float64, last :: Float64)
    last = prevfloat(last)
    firstExp = SubString(bitstring(first), 2:12)
    lastExp = SubString(bitstring(last), 2:12)
    # jeżeli eksponenty nie są równe, nie będziemy mieli równego rozmieszczenia
    if (firstExp != lastExp)
        # zero jako błąd
        return 0.0
    end

    return ((2.0^(parse(Int, firstExp, base = 2) - 1023))*2.0^(-52))
end

println(numberSpread(0.5, 1.0))
println(numberSpread(1.0, 2.0))
println(numberSpread(2.0, 4.0))