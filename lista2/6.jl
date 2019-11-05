function recur(n, c, x0)
    if (n == 0)
        return x0
    end

    return recur(n - 1, c, x0)^2 + c
end

println(recur(40, -2, 1))
println(recur(40, -2, 2))
println(recur(40, -2, 1.99999999999999))
println(recur(40, -1, 1))
println(recur(40, -1, -1))
println(recur(40, -1, 0.75))
println(recur(40, -1, 0.25))