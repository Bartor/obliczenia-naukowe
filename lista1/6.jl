function f(number)
    sqrt(number^2 + 1) - 1
end

function g(number)
    number^2/(sqrt(number^2 + 1) + 1)
end

foreach(x -> println(f(8.0^(-x)), " ", g(8.0^(-x))), 1:10)