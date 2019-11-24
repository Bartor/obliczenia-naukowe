include("./module.jl")
using .Solvers

function f(x)
    exp(x) - 3x
end

delta = 10^(-4)
epsilon = 10^(-4)

println(bisect(f, 0.0, 1.0, delta, epsilon))
println(bisect(f, 1.0, 2.0, delta, epsilon))