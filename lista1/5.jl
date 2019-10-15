# Autor: Bartosz Rajczyk

"Wektory jako liczby zmiennoprzecinkowe podwójnej precyzji"
x64 = [2.718281828, -3.141592654, 1.414213562, 0.5772156649, 0.3010299957]
y64 = [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

"Wektory jako liczby zmiennoprzecinkowe pojedynczej precyzji"
x32 = map(element -> convert(Float32, element), x64)
y32 = map(element -> convert(Float32, element), y64)

"Iloczyn skalarny dwóch wektorów \"w przód\""
function a(x, y)
    s = 0.0
    for i in (1:length(x))
        s += x[i] * y[i]
    end
    s
end

"Iloczyn skalarny dwóch wektorów \"w tył\""
function b(x, y)
    s = 0.0
    for i in reverse(1:length(x))
        s += x[i] * y[i]
    end
    s
end

"Iloczyn skalarny dwóch wektorów od dodając w kolejności od najmniejszych ujemnych/największych dodatnich kolejne sumy częściowe"
function c(x, y)
    partials = map(x -> x[1]*x[2], zip(x, y))
    foldl(+, sort(filter(x -> x > 0, partials), rev=true)) + foldl(+, sort(filter(x -> x <= 0, partials)))
end

"Iloczyn skalarny dwóch wektorów od dodając w kolejności od najmniejszych dodatnich/największych ujemnych kolejne sumy częściowe"
function d(x, y)
    partials = map(x -> x[1]*x[2], zip(x, y))
    foldl(+, sort(filter(x -> x > 0, partials))) + foldl(+, sort(filter(x -> x <= 0, partials), rev=true))
end

fs = [a, b, c, d]
println("Float64")
foreach(func -> println(func(x64, y64)), fs)
println("Float32")
foreach(func -> println(func(x32, y32)), fs)