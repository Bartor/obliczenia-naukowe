x = [2.718281828, -3.141592654, 1.414213562, 0.5772156649, 0.3010299957]
y = [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

function a()
    s = 0.0
    for i in (1:length(x))
        s += x[i] * y[i]
    end
    s
end

function b()
    s = 0.0
    for i in reverse(1:length(x))
        s += x[i] * y[i]
    end
    s
end

function c()
    partials = map(x -> x[1]*x[2], zip(x, y))
    foldl(+, sort(filter(x -> x > 0, partials), rev=true)) + foldl(+, sort(filter(x -> x <= 0, partials)))
end

function d()
    partials = map(x -> x[1]*x[2], zip(x, y))
    foldl(+, sort(filter(x -> x > 0, partials))) + foldl(+, sort(filter(x -> x <= 0, partials), rev=true))
end

println(a())
println(b())
println(c())
println(d())