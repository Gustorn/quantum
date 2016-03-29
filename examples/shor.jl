include("../src/qspice.jl")

using QSpice.Gates, QSpice.State, QSpice.BitOps

function pow2between(n1, n2)
    for i = n1:n2
        if (i & (i - 1)) == 0
            return i
        end
    end
    error("Cannot find appropriate power of 2 in the given interval")
end

# Two random integers have a ~60% chance of being coprimes, so just resamping
# until we find one is a viable algorithm
function randcoprime(n)
    cprm = rand(3:n - 1)
    while gcd(cprm, n) != 1
        cprm = rand(3:n - 1)
    end
    return cprm
end

function approxrational(n::Real, limit::Int)
    m = [1 0; 0 1]

    current = n
    cint = floor(Int, current)

    while m[2,1] * cint + m[2, 2] <= limit
        m = m * [cint 1; 1 0]
        if current == float(cint)
            break
        end
        current = 1.0 / (current - cint)
        cint = floor(Int, current)
    end

    e1 = n - m[1,1] / m[2,1]

    cint = floor(Int, (limit - m[2,2]) / m[2,1])
    n2 = m[1,1] * cint + m[1,2]
    d2 = m[2,1] * cint + m[2,2]
    e2 = n - n2 / d2

    if abs(e1) <= abs(e2)
        return m[1,1], m[2,1]
    else
        return n2, d2
    end
end

function factorize(n)
    if isprime(n)
        return [n]
    end

    if iseven(n)
        return [2, factorize(n / 2)...]
    end

    q = pow2between(n * n, 2 * n * n - 1)
    x = randcoprime(n)

    factor = 1
    while factor == 1
        m = findperiod(n, q, x)
        r = approxrational(m / q, q)[2]
        factor = max(gcd(x^div(r, 2) + 1, n) % n,
                     gcd(x^div(r, 2) - 1, n) % n)
    end

    return [factor, div(n, factor)]
end

function findperiod(n, q, x)
    # First initialize register 1 with an equally weighed superposition of integers
    # from 0 to q - 1
    reg1 = fromstates([QUBIT0 for _ in 1:floor(Int, log2(q))]...)
    for i = 1:reg1.bits
        reg1 = hadamard(reg1, i)
    end

    # This part is not yet done with actual quantum gates/functions
    bits = ceil(Int, log2(n))
    reg2 = QuantumState([0.0 + 0im for _ in 1:2^bits], bits)
    for a = 0:q - 1
        reg2[powermod(x, a, n) + 1] += reg1[a + 1]
    end
    n = sqrt(sumabs2(reg2.vector))
    reg2.vector = reg2.vector .* (1.0 / n)

    reg = Gates.join(reg1, reg2)
    post, k = partialmeasure(reg, reg1.bits + 1:reg.bits...)

    collapsed, m = measure(qft(post))
    return todecimal(m)
end
