function a()
    current = one(Float64)
    while nextfloat(current) * (one(Float64)/nextfloat(current)) == one(Float64)
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = a()
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))

function b()
    current = nextfloat(zero(Float64))
    while nextfloat(current) * (one(Float64)/nextfloat(current)) == one(Float64)
        current = nextfloat(current)
    end
    nextfloat(current)
end

result = b()
println("for ", result, " x*(1/x) = ", result*(one(Float64)/result))
