# Autor: Bartosz Rajczyk

module LinearSolvers
export matrixFromInput, vectorFromInput, calculateRightSide, gauss!, gaussWithPivots!, gaussLU!, gaussWithPivotsLU!, solveGauss, solveGaussWithPivots, solveLU, solveLUWithPivots

using SparseArrays

function matrixFromInput(file::IOStream)
    words = split(readline(file))
    matrixSize = parse(Int, words[1])
    blockSize = parse(Int, words[2])

    columns = []
    rows = []
    values = []

    for (lineNumber, line) in enumerate(Iterators.drop(eachline(file), 1))
        words = split(line)
        push!(columns, parse(Int, words[1]))
        push!(rows, parse(Int, words[2]))
        push!(values, parse(Float64, words[3]))
    end

    M = sparse(columns, rows, values)
    return (M, matrixSize, blockSize)
end

function vectorFromInput(file::IOStream)
    words = split(readline(file))
    size = parse(Int, words[1])

    vector = Array{Float64}(undef, size)

    for (lineNumber, line) in enumerate(Iterators.drop(eachline(file), 1))
        vector[lineNumber] = parse(Float64, line)
    end
        
    return vector
end

function calculateRightSide(M::SparseMatrixCSC{Float64,Int}, size::Int, blockSize::Int)
    solution = zeros(Float64, size)
    for i in 1:size
        startC = convert(Int, max(i - (2 + blockSize - i % blockSize), 1))
        endC = convert(Int, min(i + blockSize, size))

        for j in startC:endC
            solution[i] += M[i, j]
        end
    end
        
    return solution
end


# In-place gauss elimination with right-side vector
function gauss!(M! :: SparseMatrixCSC{Float64,Int}, b! :: Vector{Float64}, size :: Int, blockSize :: Int)
    for i in 1:size-1
        for j in i+1:min(size, i + 1 + blockSize)
            if eps(Float64) > abs(M![i, i])
                error("Diagol coefficient smaller than machine epsilon")
            end

            z = M![j, i] / M![j, j]
            M![j, i] = 0.0

            for k in i+1:min(size, i + blockSize + i % blockSize + 2)
                M![j, k] = M![j, k] - z * M![i, k]
            end
            b![j] = b![j] - z * b![i]
        end
    end
end

# In-place gauss elimination with choice
function gaussWithPivots!(M! :: SparseMatrixCSC{Float64,Int}, b! :: Vector{Float64}, size :: Int, blockSize :: Int) :: Vector{Int}
    pivots = collect(1:size)

    for k in 1:size-1
        lastColumn = 0
        lastRow = 0

        for i in k:min(size, k + blockSize + 1)
            if abs(M![pivots[i], k]) > lastColumn
                lastColumn = abs(M![pivots[i], k])
                lastRow = i
            end
        end

        pivots[lastRow], pivots[k] = pivots[k], pivots[lastRow]

        for i in k+1:min(size, k + blockSize + 1)
            z = M![pivots[i], k] / M![pivots[k], k]
            M![pivots[i], k] = 0.0

            for j in k+1:min(size, k + blockSize + k % blockSize + 2)
                M![pivots[i], j] = M![pivots[i], j] - z * M![pivots[k], j]
            end
            b![pivots[i]] = b![pivots[i]] - z * b![pivots[k]]
        end
    end

    return pivots
end

# where M and b are from gauss
function solveGauss(M :: SparseMatrixCSC{Float64,Int}, b :: Vector{Float64}, size :: Int, blockSize :: Int) :: Vector{Float64}
    result = zeros(Float64, size)

    for i in size:-1:1
        currentSum = 0
        for j in i+1:min(size, i + blockSize)
            currentSum += M[i, j] * result[j]
        end

        result[i] = (b[i] - currentSum) / M[i, i]
    end

    return result
end

# Where M and b are from gauss with pivots
function solveGaussWithPivots(M :: SparseMatrixCSC{Float64, Int}, b :: Vector{Float64}, size :: Int, blockSize :: Int, pivots :: Vector{Int}) :: Vector{Float64}
    result = zeros(Float64, size)
    
    for k in 1:size-1
        for i in k+1:min(size, k + 2 * blockSize)
            b[pivots[i]] = b[pivots[i]] - M[pivots[i], k] * b[pivots[k]]
        end
    end

    for i in size:-1:1
        currentSum = 0
        for j in i+1:min(size, i + 2 * blockSize)
            currentSum += M[pivots[i], j] * result[j]
        end
        result[i] = (b[pivots[i]] - currentSum) / M[pivots[i], i]
    end

    return result
end

# Where L! is original matrix, U! is zero-matrix of the same size
function gaussLU!(L! :: SparseMatrixCSC{Float64, Int}, U! :: SparseMatrixCSC{Float64, Int}, size :: Int, blockSize :: Int)
    for k in 1:size-1
        L![k, k] = 1
        for i in k+1:min(size, k + blockSize + 1)
            z = U![i, k] / U![k, k]
            L![i, k] = z
            U![i, k] = 0
            for j in k+1:min(size, k + 2 * blockSize)
                U![i, j] = U![i, j] - z *U![k, j]
            end
        end
    end
    L![size, size] = 1
end

function gaussWithPivotsLU!(U! :: SparseMatrixCSC{Float64, Int}, L! :: SparseMatrixCSC{Float64, Int}, size :: Int, blockSize :: Int) :: Vector{Int}
    pivots = collect(1:size)

    for k in 1:size-1
        maximumColumnValue = 0
        maximumIndex = 0

        for i in k:min(size, k + blockSize + 1)
            if abs(U![pivots[i], k]) > maximumColumnValue
                maximumColumnValue = abs(U![pivots[i], k])
                maximumIndex = i
            end
        end

        pivots[maximumIndex], pivots[k] = pivots[k], pivots[maximumIndex]

        for i in k+1:min(size, k + blockSize + 1)
            z = U![pivots[i], k] / U![pivots[k], k]

            L![pivots[i], k] = z
            U![pivots[i], k] = 0

            for j in k+1:min(size, k + 2 * blockSize)
                U![pivots[i], j] = U![pivots[i], j] - z * U![pivots[k], j]
            end
        end
    end

    return pivots
end

function solveLU(L :: SparseMatrixCSC{Float64,Int}, U :: SparseMatrixCSC{Float64,Int}, b :: Vector{Float64}, size :: Int, blockSize :: Int)
    result = zeros(Float64, size)

    for k in 1:size-1
        for i in k+1:min(size, k + blockSize + 1)
            b[i] = b[i] - L[i, k] * b[k]
        end
    end

    for i in size:-1:1
        currentSum = 0
        for j in i+1:min(size, i + blockSize)
            currentSum += U[i, j] * result[j]
        end
        result[i] = (b[i] - currentSum) / U[i, i]
    end

    return result
end

function solveLUWithPivots(L :: SparseMatrixCSC{Float64,Int}, U :: SparseMatrixCSC{Float64,Int}, b :: Vector{Float64}, size :: Int, blockSize :: Int, pivots :: Vector{Int})
    result = zeros(Float64, size)

    for k in 1:size-1
        for i in k+1:min(size, 2 * blockSize + k + 5)
            b[pivots[i]] = b[pivots[i]] - L[pivots[i], k] * b[pivots[k]]
        end
    end

    for i in size:-1:1
        currentSum = 0
        for j in i+1:min(size, i + 2 * blockSize)
            currentSum += U[pivots[i], j] * result[j]
        end
        result[i] = (b[pivots[i]] - currentSum) / U[pivots[i], i]
    end

    return result
end

# module end
end