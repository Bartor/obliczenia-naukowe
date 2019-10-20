# Autor: Bartosz Rajczyk

function deriviative(fun, x, h)
    (fun(x + h) - fun(x))/h
end

function fun(x)
    sin(x) + cos(3*x)
end

function funPrim(x)
    cos(x) - 3*sin(3*x)
end

result = funPrim(1.0)
foreach(h -> println("\$2^{-$h}\$", " & ", 1.0 + 2.0^(-h) ," & ", deriviative(fun, 1.0, 2.0^(-h)), " & ", abs(result - deriviative(fun, 1.0, 2.0^(-h))), " \\\\\n\\hline"), 0:54)
println(result)