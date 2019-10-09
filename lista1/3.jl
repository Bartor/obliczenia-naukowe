function floatRepresentation(first, last)
    println(bitstring(first))
    println(bitstring(last))
end

floatRepresentation(0.5, prevfloat(1.0))