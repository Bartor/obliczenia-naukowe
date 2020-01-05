include("./matrixgen.jl")
include("./module.jl")
using .LinearSolvers
using .matrixgen
using SparseArrays
# using PyPlot

REP_NUMER = 25
MAX_SIZE = 15000
JUMP_SIZE = 400
BLOCK_SIZE = 4

function challenge(M, L, b, size, blockSize)
    p = gaussWithPivotsLU!(M, L, size, blockSize)
    solveLUWithPivots(M, L, b, size, blockSize, p)
end

function benchmark()
    for size in JUMP_SIZE:JUMP_SIZE:MAX_SIZE
        totalTime = 0
        totalMemory = 0
        for reps in 1:REP_NUMER
            M = blockmat(size, BLOCK_SIZE, 1.0)
            b = calculateRightSide(M, size, BLOCK_SIZE)
            L = SparseArrays.spzeros(size, size)
            (_, time, memory) = @timed challenge(M, L, b, size, BLOCK_SIZE)
            totalTime += time
            totalMemory += memory
        end
        println(size, "; ", replace(string(totalTime / REP_NUMER), "." => ","), "; ", replace(string(totalMemory / REP_NUMER), "." => ","))
    end
end

# M = blockmat(50000, BLOCK_SIZE, 1.0)
# b = calculateRightSide(M, 50000, BLOCK_SIZE)
# L = SparseArrays.spzeros(50000, 50000)
# (_, time, memory) = @timed challenge(M, L, b, 50000, BLOCK_SIZE)

# println(time)
# println(memory)

benchmark()

# size = 16
# blockSize = 4
# matrixFile = open("data/16_1_1/A.txt")
# vectorFile = open("data/16_1_1/b.txt")
# # M = matrixFromInput(matrixFile)
# # b = vectorFromInput(vectorFile)
# M = blockmat(size, blockSize, 1.0)
# b = calculateRightSide(M, size, blockSize)
# L = SparseArrays.spzeros(size, size)
# p = gaussWithPivotsLU!(M, L, size, blockSize)
# solution = solveLUWithPivots(M, L, b, size, blockSize, p)
# println(solution)