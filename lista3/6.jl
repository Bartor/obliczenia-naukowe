include("./module.jl")
using .Solvers

function f1(x)
    exp(1-x) - 1
end

function f1Prim(x)
    -exp(1-x)
end

function f2(x)
    x*exp(-x)
end

function f2Prim(x)
    -exp(-x)*(x-1)
end

delta = 10^(-5)
epsilon = 10^(-5)

println(bisect(f1, 0.0, 2.0, delta, epsilon))
println(bisect(f2, -1.0, 1.0, delta, epsilon))

println(newton(f1, f1Prim, 0.0, delta, epsilon, 100))
println(newton(f2, f2Prim, -1.0, delta, epsilon, 100))

println(secant(f1, 0.0, 0.1, delta, epsilon, 100))
println(secant(f2, -1.0, -0.9, delta, epsilon, 100))

println(newton(f1, f1Prim, 100.0, delta, epsilon, 100))
println(newton(f2, f2Prim, 100.0, delta, epsilon, 100))
println(newton(f2, f2Prim, 1.0, delta, epsilon, 100))