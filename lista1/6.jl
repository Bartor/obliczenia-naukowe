# Autor: Bartosz Rajczyk

function f(number)
    sqrt((number^2) + 1) - 1
end

function g(number)
    (number^2)/(sqrt((number^2) + 1) + 1)
end

# Obliczanie obu funkcji od 8^(-1) do 8^(-180) - granica wyznaczona eksperymentalnie
foreach(x -> println(x, " & ", f(8.0^(-x)), " & ", g(8.0^(-x)), "\\\\\n\\hline"), 1:9)
foreach(x -> println(20*x, " & ", f(8.0^(-20*x)), " & ", g(8.0^(-20*x)), "\\\\\n\\hline"), 1:9)