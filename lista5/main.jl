include("./matrixgen.jl")
include("./module.jl")
using .LinearSolvers
using .matrixgen
using SparseArrays
# using PyPlot

REP_NUMER = 10
MAX_SIZE = 15000
JUMP_SIZE = 400
BLOCK_SIZE = 4

function challange(M, L, b, size, blockSize)
    gauss!(M, b, size, blockSize)
    solveGauss(M, b, size, blockSize)
end

function benchmark()
    for size in JUMP_SIZE:JUMP_SIZE:MAX_SIZE
        totalTime = 0
        totalMemory = 0
        for reps in 1:REP_NUMER
            M = blockmat(size, BLOCK_SIZE, 1.0)
            b = calculateRightSide(M, size, BLOCK_SIZE)
            L = SparseArrays.spzeros(size, size)
            (_, time, memory) = @timed challange(M, L, b, size, BLOCK_SIZE)
            totalTime += time
            totalMemory += memory
        end
        println(size, ", ", totalTime / REP_NUMER, ", ", totalMemory / REP_NUMER)
    end
end

benchmark()

#=

matrixFile = open("data/16_1_1/A.txt")
vectorFile = open("data/16_1_1/b.txt")
(M, size, blockSize) = matrixFromInput(matrixFile)
b = vectorFromInput(vectorFile)
pivots = matrixLUWithPivots!(M, size, blockSize)
solution = solveLUWithPivots(M, b, size, blockSize, pivots)
println(solution)

=#