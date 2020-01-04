# Autor: Bartosz Rajczyk

module LinearSolvers
export matrixFromInput, vectorFromInput, gaussianElimination, gaussianWithPivots, matrixLU!, matrixLUWithPivots!, solveLU, solveLUWithPivots

using SparseArrays

function matrixFromInput(file::IOStream)
    words = split(readline(file))
    matrixSize = parse(Int, words[1])
    blockSize = parse(Int, words[2])
    blockNumber = convert(Int, matrixSize / blockSize)
    arraySize = blockNumber * blockSize * blockSize + 3 * (blockNumber - 1) * blockSize - 1

    columns = Array{Int}(undef, arraySize)
    rows = Array{Int}(undef, arraySize)
    values = Array{Float64}(undef, arraySize)

    for (lineNumber, line) in enumerate(Iterators.drop(eachline(file), 1))
        words = split(line)
        columns[lineNumber] = parse(Int, words[1])
        rows[lineNumber] = parse(Int, words[2])
        values[lineNumber] = parse(Float64, words[3])
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
        startC = convert(Int, max(blockSize * floor((i - 1) / blockSize) - 1, 1))
        endC = convert(Int, min(blockSize + blockSize * floor((i - 1) / blockSize), size))

        for j in startC:endC
            solution[i] += M[j, i]
        end

        if (i + blockSize > size)
            solution[i] += M[i + blockSize, i]
        end
    end
        
    return solution
end

function gaussianElimination(M::SparseMatrixCSC{Float64,Int}, b::Vector{Float64}, size::Int, blockSize::Int)
    for i in 1:size - 1
        lastColumn = convert(Int, min(i + blockSize, size))
        lastRow = convert(Int, min(blockSize + blockSize * floor((i + 1) / blockSize), size))

        for j in i + 1:lastRow
            if eps(Float64) > abs(M[i, i])
                error("Diagol coefficient smaller than machine epsilon")
            end

            z = M[i, j] / M[i, i]
            M[i, j] = 0
            for k in i + 1:lastColumn
                M[k, j] = M[k, j] - z * M[k, i]
            end
            b[j] = b[j] - z * b[i]
        end
    end

    result = Array{Float64}(undef, size)
    for i in size:-1:1
        sum = 0
        last = min(size, i + 1)
        for j in i + 1:last
            sum += M[j, i] * result[j]
        end
        result[j] = (b[i] - sum) / M[i, i]
    end

    return result
end

function gaussianWithPivots(M::SparseMatrixCSC{Float64,Int}, b::Vector{Float64}, size::Int, blockSize::Int)
    pivots = collect(1:size)

    for k in 1:size - 1
        lastColumn = convert(Int, min(2 * blockSize + blockSize * floor((k + 1) / blockSize), size))
        lastRow = convert(Int, min(blockSize + blockSize * floor((k + 1) / blockSize), size))
        for i in k + 1:lastRow
            maxRow = k
            max = abs(M[k, pivots[k]])
            for x in i:lastRow
                if (abs(M[k, pivots[x]]) > max)
                    maxRow = x
                    max = abs(M[k, pivots[x]])
                end
            end

            if eps(Float64) > abs(max)
                error("Matrix is singular")
            end

            pivots[k], pivots[maxRow] = pivots[maxRow], pivots[k]
            z = M[k, pivots[i]] / M[k, pivots[k]]
            M[k, pivots[i]] = 0

            for j in k + 1:lastColumn
                M[j, pivots[i]] = M[j, pivots[i]] - z * M[j, pivots[k]]
            end
            b[pivots[i]] = b[pivots[i]] - z * b[pivots[k]]
        end
    end

    result = Array{Float64}(undef, size)
    for i in size:-1:1
        lastColumn = convert(Int, min(2 * blockSize + blockSize * floor((pivots[i] + 1) / blockSize), size))
        sum = 0
        for j in i + 1:lastColumn
            sum += M[j, pivots[i]] * result[j]
        end
        result[i] = (b[pivots[i]] - sum) / M[i, pivots[i]]
    end

    return result
end

function matrixLU!(M::SparseMatrixCSC{Float64,Int}, size::Int, blockSize::Int)
    for k in 1:size - 1
        lastRow = convert(Int, min(blockSize + blockSize * floor((k + 1) / blockSize), size))
        lastColumn = convert(Int, min(blockSize + k, size))

        for i in k + 1:lastRow
            if abs(M[k, k]) < eps(Float64)
                error("Diagol coefficient smaller than machine epsilon")
            end
            z = M[k, i] / M[k, k]
            M[k, i] = z
            for j in k + 1:lastColumn
                M[j, i] = M[j, i] - z * M[j, k]
            end
        end
    end
end

function matrixLUWithPivots!(M::SparseMatrixCSC{Float64,Int}, size::Int, blockSize::Int)
    pivots = collect(1:size)
        
    for k in 1:size - 1
        lastRow = convert(Int, min(blockSize + blockSize * floor((k + 1) / blockSize), size))
        lastColumn = convert(Int, min(2 * blockSize + blockSize * floor((k + 1) / blockSize), size))

        for i in k + 1:lastRow
            maxRow = k
            max = abs(M[k, pivots[k]])
            for x in i:lastRow
                if abs(M[k, pivots[x]]) > max
                    maxRow = x
                    max = abs(M[k, pivots[x]])
                end
            end
            if eps(Float64) > abs(max)
                error("Matrix is singular")
            end
            pivots[k], pivots[maxRow] = pivots[maxRow], pivots[k]
            z = M[k, pivots[i]] / M[k, pivots[k]]
            M[k, pivots[i]] = z
            for j in k + 1:lastColumn
                M[j, pivots[i]] = M[j, pivots[i]] - z * M[j, pivots[k]]
            end
        end
    end

    return pivots
end

function solveLU(M::SparseMatrixCSC{Float64,Int}, b::Vector{Float64}, size::Int, blockSize::Int) 
    z = Array{Float64}(undef, size)
    for i in 1:size
        sum = 0
        column = convert(Int, max(blockSize * floor((i = 1) / blockSize) - 1, 1))
        for j in column:i - 1
            sum += M[j, i] * z[j]
        end
        z[i] = b[i] - sum
    end

    result = Array{Float64}(undef, size)
    for i in size:-1:1
        sum = 0
        lastColumn = min(size, i + 1)
        for j in i + 1:lastColumn
            sum += M[j, i] * result[j]
        end
        result[i] = (z[i] - sum) / M[i, i]
    end

    return result
end

function solveLUWithPivots(M::SparseMatrixCSC{Float64,Int}, b::Vector{Float64}, size::Int, blockSize::Int, pivots :: Vector{Int})
    z = Array{Float64}(undef, size)

    for i in 1:size
        sum = 0
        startColumn = convert(Int, max(blockSize * floor((i - 1) / blockSize) - 1, 1))
        for j in startColumn:i-1
            sum += M[j, pivots[i]] * z[j]
        end
        z[i] = b[pivots[i]] - sum
    end

    result = Array{Float64}(undef, size)
    for i in size:-1:1
        sum = 0
        lastColumn = convert(Int, min(2 * blockSize + blockSize * floor((pivots[i] + 1) / blockSize), size))
        for j in i+1:lastColumn
            sum += M[j, pivots[i]] * result[j]
        end
        result[i] = (z[i] - sum) / M[i, pivots[i]]
    end

    return result
end

# module end
end