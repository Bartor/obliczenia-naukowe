include("./module.jl")
using .LinearSolvers

matrixFile = open("data/10000_1_1/A.txt")
vectorFile = open("data/10000_1_1/b.txt")
(M, size, blockSize) = matrixFromInput(matrixFile)
b = vectorFromInput(vectorFile)
pivots = matrixLUWithPivots!(M, size, blockSize)
solution = solveLUWithPivots(M, b, size, blockSize, pivots)
println(solution)