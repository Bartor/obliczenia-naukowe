# Autor: Bartosz Rajczyk

include("./module.jl")
using .Interpolation

e(x) = exp(x)
x2sin(x) = x^2 * sin(x)

drawFunction(e, 0.0, 1.0, 5)
drawFunction(e, 0.0, 1.0, 10)
drawFunction(e, 0.0, 1.0, 15)

drawFunction(x2sin, -1.0, 1.0, 5)
drawFunction(x2sin, -1.0, 1.0, 10)
drawFunction(x2sin, -1.0, 1.0, 15)