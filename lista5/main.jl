include("./matrixgen.jl")
include("./module.jl")
using .LinearSolvers
using .matrixgen
# using PyPlot

REP_NUMER = 100
MAX_SIZE = 40000
JUMP_SIZE = 400

function challange(M, b, size, blockSize)
    pivots = matrixLUWithPivots!(M, size, blockSize)
    solveLUWithPivots(M, b, size, blockSize, pivots)
end

function benchmark()
    for size in 4:JUMP_SIZE:MAX_SIZE
        totalTime = 0
        totalMemory = 0
        for reps in 1:REP_NUMER
            # old way of testing
            # blockmat(size, 4, 1.0, "test")
            # file = open("test")
            # (M, size, blockSize) = matrixFromInput(file)
            # close(file)
            M = blockmat(size, 4, 1.0)
            b = calculateRightSide(M, size, 4)
            (_, time, memory) = @timed challange(M, b, size, 4)
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