export simps, trapz

function simps(xs, ys)
    N = length(xs) - 1
    h = diff(xs)
    I = 0.0
    for i in 2:2:N
        hph = h[i] + h[i-1]
        I += ys[i]   * (h[i]^3 + h[i-1]^3 + 3*h[i]*h[i-1]*hph)/(6*h[i]*h[i-1])
        I += ys[i-1] * (2*h[i-1]^3 - h[i]^3 + 3*h[i]*h[i-1]^2)/(6*h[i-1]*hph)
        I += ys[i+1] * (2*h[i]^3 - h[i-1]^3 + 3*h[i-1]*h[i]^2)/(6*h[i]*hph)
    end

    if (N + 1) % 2 == 0
        I += ys[N+1]   * (2*h[N]^2 + 3*h[N-1] * h[N])/(6*(h[N-1] + h[N]))
        I += ys[N]     * (h[N]^2 + 3*h[N]* h[N-1])/(6*h[N-1])
        I -= ys[N-1]   * (h[N]^3)/(6*h[N-1]*(h[N-1] + h[N]))
    end
    return I
end

function trapz(xs, ys)
    # number of intervals
    N = length(xs)-1
    length(ys) == length(xs) || throw(ArgumentError("argh!!!"))
    I = (ys[2]+ys[1])*(xs[2]-xs[1])
    for k = 2:N
        I += (ys[k+1]+ys[k])*(xs[k+1]-xs[k])
    end
    return I/2
end