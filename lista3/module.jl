module Solvers

export bisect, newton, secant

function bisect(f :: Function, a :: Float64, b :: Float64, delta :: Float64, epsilon :: Float64)
    if (f(a)*f(b) >= 0)
        return (0, 0, 0, 1)
    end

    it = 0
    
    e = b - a; u = f(a); v = f(b)
    while true
        it = it + 1
        e = e/2
        c = a + e
        w = f(c)

        if (abs(e) < delta || abs(w) < epsilon)
            return (e, w, it, 0)
        end

        if (w*u < 0)
            b = c;
            v = w
        else
            a = c
            u = w
        end
    end
end

function newton(f :: Function, fPrim :: Function, x0 :: Float64, delta :: Float64, epsilon :: Float64, maxIt :: Int)
    v = f(x0)
    if (abs(v) < epsilon)
        return (x0, v, 0, 0)
    end

    for i in 1:maxIt
        x1 = x0 - v/fPrim(x0)
        v = f(x1)
        if (abs(x1 - x0) < delta || abs(v) < epsilon)
            return (x1, v, i, 0)
        end
        x0 = x1
    end
    return (0, 0, 0, 1)
end

function secant(f :: Function, a :: Float64, b :: Float64, delta :: Float64, epsilon :: Float64, maxIt :: Int)
    fa = f(a)
    fb = f(b)

    for i in 1:maxIt
        if (abs(fa) > abs(fb))
            a, b = b, a
            fa, fb = fb, fa
        end
        s = (b - a)/(fb - fa)
        b = a
        fb = fa
        a = a - fa*s
        fa = f(a)
        if (abs(b - a) < delta || abs(fa) < epsilon)
            return (a, fa, i, 0)
        end
    end
    return (0, 0, 0, 1)
end

end