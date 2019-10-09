types = [Float16, Float32, Float64]

function machEps(type)
    number = one(type)
    while one(number) + number / 2 != one(number)
        number /= 2
    end
    number
end

println("Machine Epsilon")
foreach(type -> println("f: ", machEps(type), " eps: ", eps(type)), types)

function eta(type)
    number = one(type)
    while number / 2 != 0
        number /= 2
    end
    number
end

println("Eta")
foreach(type -> println("f: ", eta(type), " nextfloat: ", nextfloat(zero(type))), types)

function maxNumber(type)
    number = one(type)
    while !isinf(number * 2)
        number *= 2
    end
    number
end

println("Max")
foreach(type -> println("f: ", maxNumber(type), " floatmax: ", floatmax(type)), types)