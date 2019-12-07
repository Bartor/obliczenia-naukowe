# Autor: Bartosz Rajczyk

include("./module.jl")
using .Interpolation

oneOverOnePlusXSquared(x) = 1/(1+x^2)

drawFunction(abs, -1.0, 1.0, 5)
drawFunction(abs, -1.0, 1.0, 10)
drawFunction(abs, -1.0, 1.0, 15)

drawFunction(oneOverOnePlusXSquared, -5.0, 5.0, 5)
drawFunction(oneOverOnePlusXSquared, -5.0, 5.0, 10)
drawFunction(oneOverOnePlusXSquared, -5.0, 5.0, 15)