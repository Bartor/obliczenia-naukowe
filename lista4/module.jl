# Autor: Bartosz Rajczyk

module Interpolation
using PyPlot

    export quotients, newtonValue, naturalForm, drawFunction

    function quotients(xs :: Vector{Float64}, fxs :: Vector{Float64})
        len = length(fxs)
        if len != length(xs)
            error("Lenghts are not equal")
        end

        res = [value for value in fxs]
        for i in 1:len
            for j in len:-1:i+1
                res[j] = (res[j] - res[j-1]) / (xs[j] - xs[j - i])
            end
        end

        return res
    end

    function newtonValue(xs :: Vector{Float64}, qs :: Vector{Float64}, t :: Float64)
        len = length(xs)
        if len != length(qs)
            error("Lengths are not equal")
        end

       nth = qs[len]
       for i in len-1:-1:1
            nth = qs[i] + (t - xs[i]) * nth
       end

       return nth
    end

    function naturalForm(xs :: Vector{Float64}, fxs :: Vector{Float64})
        len = length(xs)
        if len != length(fxs)
            error("Lenghts are not equal")
        end

        res = [value for value in fxs]
        for i in len-1:-1:1
            res[i] = fxs[i] - res[i+1] * xs[i]
            for j = i+1:len-1
                res[j] = res[j]-res[j+1]*xs[i]
            end
        end
        return res
    end

    function drawFunction(f, a :: Float64, b :: Float64, n :: Int)
        if a > b
            a, b = b, a
        end

        steps = n + 1

        xs = Vector{Float64}(undef, steps)
        ys = Vector{Float64}(undef, steps)

        currentDelta = zero(Float64)
        for i in 1:steps
            xs[i] = a + currentDelta
            ys[i] = f(xs[i])
            currentDelta += (b - a) / (steps - 1)
        end

        density = 10
        steps = (n + 1) * density

        qs = quotients(xs, ys)
        interpolationXs = Vector{Float64}(undef, steps)
        interpolationVals = Vector{Float64}(undef, steps)
        realVals = Vector{Float64}(undef, steps)

        currentDelta = zero(Float64)
        for i in 1:steps
            interpolationXs[i] = a + currentDelta
            interpolationVals[i] = newtonValue(xs, qs, a + currentDelta)
            realVals[i] = f(a + currentDelta)
            currentDelta += (b - a) / (steps - 1)
        end

        clf()
        plot(interpolationXs, interpolationVals, label="interpolated", linewidth=2.0, alpha=0.5, color="#0070ff")
        plot(interpolationXs, realVals, label="actual", linewidth=2.0, alpha=0.5, color="#ff7000")
        grid(true)
        legend(title="Interpolation")
        savefig(string("plotted/plot", f, "-", n, ".png"))
    end
end
