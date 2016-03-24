module Gates

using QSpice.BitOps
using QSpice.State

export qidentity, hadamard, not, cnot, swap, sqrtswap,
       phaseshift, paulix, pauliy, pauliz,
       ccnot, cswap, unitary, choose1,
       measure, partialmeasure, probe, superposition, qft

function qidentity(state::QuantumState)
    return copy(state)
end

function superposition(state1::QuantumState, state2::QuantumState, ss::QuantumState...)
    return fromstates(state1, state2, ss...)
end

function hadamard(state::QuantumState, bit::Int)
    newstate = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit
    invsqrt2 = 1.0 / sqrt(2)

    for i = 1:length(state)
        basis  = i - 1
        mapped0 = clearbit(basis, bit)
        mapped1 = setbit(basis, bit)

        if iszero(basis, bit)
            newstate[i] = invsqrt2 * (state[mapped0 + 1] + state[mapped1 + 1])
        else
            newstate[i] = invsqrt2 * (state[mapped0 + 1] - state[mapped1 + 1])
        end
    end

    return QuantumState(newstate, state.bits)
end

function not(state::QuantumState, bit::Int)
    newstate = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        mapped = flipbit(basis, bit)
        newstate[i] = state[mapped + 1]
    end

    return QuantumState(newstate, state.bits)
end

function cnot(state::QuantumState, ctrl::Int, flip::Int)
    newstate = Array{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    flip = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, ctrl)
            newstate[i] = state[i]
        else
            mapped = flipbit(basis, flip)
            newstate[i] = state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

function ccnot(state::QuantumState, ctrl1::Int, ctrl2::Int, flip::Int)
    newstate = Array{Complex{Float64}}(length(state))

    ctrl1 = state.bits - ctrl1
    ctrl2 = state.bits - ctrl2
    flip  = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if !iszero(basis, ctrl1) && !iszero(basis, ctrl2)
            mapped = flipbit(basis, flip)
            newstate[i] = state[mapped + 1]
        else
            newstate[i] = state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

function swap(state::QuantumState, x::Int, y::Int)
    newstate = Array{Complex{Float64}}(length(state))

    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        mapped = swapbit(basis, x, y)
        newstate[i] = state[mapped + 1]
    end

    return QuantumState(newstate, state.bits)
end

function cswap(state::QuantumState, ctrl::Int, x::Int, y::Int)
    newstate = Array{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, ctrl)
            newstate[i] = state[i]
        else
            mapped = swapbit(basis, x, y)
            newstate[i] = state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

function sqrtswap(state::QuantumState, x::Int, y::Int)
    newstate = Array{Complex{Float64}}(length(state))

    x = state.bits - x
    y = state.bits - y

    mulorig = 0.5 + 0.5im
    mulswap = 0.5 - 0.5im

    for i = 1:length(state)
        basis = i - 1
        mapped = swapbit(basis, x, y)
        if basis == mapped
            newstate[i] = state[i]
        else
            newstate[i] = mulorig * state[i] + mulswap * state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

function phaseshift(state::QuantumState, bit::Int, theta::Float64)
    newstate = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit
    mul = exp(im * theta)

    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, bit)
            newstate[i] = state[i]
        else
            newstate[i] = mul * state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

paulix(state::QuantumState, bit::Int) = not(state, bit)

function pauliy(state::QuantumState, bit::Int)
    newstate = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, bit)
            mapped = setbit(basis, bit)
            newstate[i] =  im * state[mapped + 1]
        else
            mapped = clearbit(basis, bit)
            newstate[i] = -im * state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

function pauliz(state::QuantumState, bit::Int)
    newstate = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, bit)
            newstate[i] =  state[i]
        else
            newstate[i] = -state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

function unitary(state::QuantumState, matrix::Array)
    if length(state.vector) != size(matrix, 1) || length(state.vector) != size(matrix, 2)
        error("The unitary matrix dimensions don't match the quantum state vector's");
    end

    return QuantumState(matrix * state.vector, state.bits)
end

function unitary(state::QuantumState, matrix::Array, x0::Int, xs::Int...)
    bits = [x0, xs...]

    if 2^length(bits) != size(matrix, 1) || 2^length(bits) != size(matrix, 2)
        error("The unitary matrix does not map enough bits")
    end

    newstate = state

    for i = 1:length(bits)
        newstate = swap(newstate, i, bits[i])

        for j = i + 1:length(bits)
            if bits[j] == i
                bits[j] = bits[i]
            end
        end
    end

    if length(state) == size(matrix, 1)
        newstate.vector = matrix * newstate.vector
    else
        gatematrix = kron(matrix, speye(2^(state.bits - length(bits))))
        newstate.vector = gatematrix * newstate.vector
    end

    sb = length(bits) + 1
    reverse!(bits)
    for i = 1:length(bits)
        newstate = swap(newstate, sb - i, bits[i])
    end

    return newstate
end

function qft(state::QuantumState)
    omega = exp(2 * pi * im / state.bits)
    matrix = [omega^((i - 1) * (j - 1)) for i = 1:length(state), j = 1:length(state)] .* (1.0 / sqrt(length(state)))
    newstate = matrix * state.vector
    return QuantumState(newstate, state.bits)
end

function qft(state::QuantumState, bit::Int, bits::Int...)
    partialsize = 2^length([bit, bits...])
    omega = exp(2 * pi * im / state.bits)
    matrix = [omega^((i - 1) * (j - 1)) for i = 1:partialsize, j = 1:partialsize] .* (1.0 / sqrt(partialsize))
    return unitary(state, matrix, bit, bits...)
end

function choose1(state::QuantumState, bits::Vector{Int}, gates::Vector{Vector{Tuple{Function, Vector{Any}}}})
    index = todecimal(bits) + 1
    chain = gates[index]

    newstate::QuantumState = chain[1][1](state, chain[1][2]...)
    for i = 2:length(chain)
        newstate = chain[i][1](newstate, chain[i][2]...)
    end
    return newstate
end

function randbasis{T<:Real}(probabilities::Vector{T})
    psum = sum(probabilities)

    target = rand() * psum
    for (i, v) in enumerate(probabilities)
        target -= v
        if (target <= 0)
            return i - 1
        end
    end

    # In case of floating point rounding error, return the element
    # with the otherwise highest probability
    return indmax(probabilities) - 1
end


function measure(state::QuantumState)
    probabilities = map(abs2, state.vector)::Vector{Float64}
    basis = randbasis(probabilities)

    # Measure the quantum state
    measurement = reverse(digits(basis, 2, state.bits))

    # Collapse the wave function
    if state.bits > 1
        newstate = fill(0.0 + 0im, length(state))
        newstate[basis + 1] = 1.0 + 0im
        posterior = QuantumState(newstate, state.bits)
    else
        posterior = EMPTY_STATE
    end

    return posterior, measurement
end

# Performs partial measurement on the given bit of the quantum state
# in the computational basis
function partialmeasure(state::QuantumState, bit::Int)
    if state.bits == 1
        return measure(state)
    end

    posterior::Vector{Complex{Float64}} = fill(0.0 + 0im, 2^(state.bits - 1))
    bit = state.bits - bit

    # First collect the coefficients for the probabilities of the bit being |0>
    prob0 = 0.0
    for i = 1:length(state)
        basis = i - 1
        if iszero(basis, bit)
            prob0 += abs2(state[i])
        end
    end

    # And measure the bit.
    measurement = rand() <= prob0 ? 0 : 1

    # Next collect the coefficients for the basis vectors of the
    # posterior states
    for i = 1:length(state)
        basis = i - 1
        if getbit(basis, bit) == measurement
            postbasis = removebit(basis, bit)
            posterior[postbasis + 1] += state[i]
        end
    end

    # Now normalize the new state
    # TODO(gustorn): see if the error from FP precision has any visible
    # effect on the posterior states
    n = sqrt(sumabs2(posterior))
    posterior = posterior .* (1.0 / n)

    return QuantumState(posterior, state.bits - 1), [measurement]
end

function partialmeasure(state::QuantumState, x0::Int, x1::Int, xs::Int...)
    # First construct the array of bits we want to measure.
    bits = [x0, x1, xs...]
    bits = unique(bits)
    nbits  = length(bits)

    if nbits == state.bits
        return measure(state)
    end

    map!(b -> state.bits - b, bits)
    sort!(bits)

    nposterior = state.bits - nbits
    probabilities::Vector{Float64} = fill(0.0, 2^nbits)
    posterior::Vector{Complex{Float64}} = fill(0.0 + 0im, 2^nposterior)

    # First build up the probabilities for the bits we're going to measure
    for i = 1:length(state)
        basis = i - 1
        measuredbasis = 0
        shiftby = 0

        for bit in bits
            measuredbasis = measuredbasis | (getbit(basis, bit) << shiftby)
            shiftby += 1
        end

        probabilities[measuredbasis + 1] += abs2(state[i])
    end

    # Measure it.
    measurement = digits(randbasis(probabilities), 2, nbits)

    # Now build up the posterior state
    for basis = 0:length(state) - 1
        # Check if the current basis vector matches the measurement at the specified bits
        basismatches = true
        for (bit, result) in zip(bits, measurement)
            if getbit(basis, bit) != result
                basismatches = false
                break
            end
        end

        # We only need to bother with the coefficients if the measurement matches the
        # partial basis vector
        if basismatches
            nremoved = 0
            postindex = basis
            for bit in bits
                postindex = removebit(postindex, bit - nremoved)
                nremoved += 1
            end

            # Add the coefficients from the old state to the new one
            posterior[postindex + 1] += state[basis + 1]
        end
    end

    # Normalize the new state
    n = sqrt(sumabs2(posterior))
    posterior = posterior .* (1 / n)

    # Match the bit order of the regular `measure` function
    reverse!(measurement)
    return QuantumState(posterior, state.bits - length(bits)), measurement
end

# Prints the current state of the qubit in terms of basis vectors, coefficients
# and probabilities. Useful for debugging and visualization. Guaranteed not to disturb
# the quantum state
function probe(state::QuantumState, name::AbstractString)
    println(name, "\n", state)
    return state
end

# The same as probe, but without the need for an explicit name
function probe(state::Vector{Int}, name::AbstractString)
    println(name, "\n", state, "\n")
    return state
end

end


