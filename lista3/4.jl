include("./module.jl")
using .Solvers

function f(x)
    sin(x) - (1/2)x^2
end

function fPrim(x)
    cos(x) - x
end

delta = (1/2)*10^(-5)
epsilon = (1/2)*10^(-5)

println(bisect(f, 1.5, 2.0, delta, epsilon))
println(newton(f, fPrim, 1.5, delta, epsilon, 100))
println(secant(f, 1.0, 2.0, delta, epsilon, 100))