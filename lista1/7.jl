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
foreach(h -> println(deriviative(fun, 1.0, 2.0^(-h)), " ", abs(result - deriviative(fun, 1.0, 2.0^(-h)))), 0:54)