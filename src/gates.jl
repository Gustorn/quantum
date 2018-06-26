module Gates

using QSpice.BitOps
using QSpice.State

export identity, hadamard, not, cnot, ccnot,
       swap, cswap, sqrtswap, phaseshift, cphase,
       paulix, pauliy, pauliz,
       unitary, choose1, superposition, joinstate, splitstate, qft,
       measure, partialmeasure, probe

function identity(state::QuantumState)
    return copy(state)
end

function superposition(state1::QuantumState, state2::QuantumState)
    return fromstates(state1, state2)
end

function joinstate(state1::QuantumState, state2::QuantumState)
    return QuantumState(vcat(state1.vector, state2.vector), state1.bits + state2.bits)
end

function splitstate(state::QuantumState, lastbit::Int)
    lpos = 2^lastbit
    return QuantumState(state[1:lpos], lastbit), QuantumState(state[lpos:end], state.bits - lastbit)
end

# Performs the following mapping:
# |0> to (|0> + |1>) / sqrt(2)
# |1> to (|0> - |1>) / sqrt(2)
function hadamard(state::QuantumState, bit::Int)
    newstate = Vector{Complex{Float64}}(length(state))
    bit = state.bits - bit
    invsqrt2 = 1.0 / sqrt(2)

    for i = 1:length(state)
        basis  = i - 1
        mapped0 = clearbit(basis, bit)
        mapped1 = setbit(basis, bit)

        if iszerobit(basis, bit)
            newstate[i] = invsqrt2 * (state[mapped0 + 1] + state[mapped1 + 1])
        else
            newstate[i] = invsqrt2 * (state[mapped0 + 1] - state[mapped1 + 1])
        end
    end

    return QuantumState(newstate, state.bits)
end

# Also called the Pauli-X gate. Performs the following mapping:
# |0> to |1>
# |1> to |0>
function not(state::QuantumState, bit::Int)
    newstate = Vector{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        mapped = flipbit(basis, bit)
        newstate[i] = state[mapped + 1]
    end

    return QuantumState(newstate, state.bits)
end

# Performs the following mapping: |x, y> to |x, x + y mod2>
function cnot(state::QuantumState, ctrl::Int, flip::Int)
    newstate = Vector{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    flip = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, ctrl)
            newstate[i] = state[i]
        else
            mapped = flipbit(basis, flip)
            newstate[i] = state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

# Performs the following mapping: |x, y, z> to |x, y, xy + z mod2>
function ccnot(state::QuantumState, ctrl1::Int, ctrl2::Int, flip::Int)
    newstate = Vector{Complex{Float64}}(length(state))

    ctrl1 = state.bits - ctrl1
    ctrl2 = state.bits - ctrl2
    flip  = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if !iszerobit(basis, ctrl1) && !iszerobit(basis, ctrl2)
            mapped = flipbit(basis, flip)
            newstate[i] = state[mapped + 1]
        else
            newstate[i] = state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

# Performs the following mapping: |x, y> to |y, x>
function swap(state::QuantumState, x::Int, y::Int)
    newstate = Vector{Complex{Float64}}(length(state))

    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        mapped = swapbit(basis, x, y)
        newstate[i] = state[mapped + 1]
    end

    return QuantumState(newstate, state.bits)
end

# Given the state |x, y, z> performs the following mapping:
# If x equals to 0
#     |x, y, z> to |x, y, z>
# If x equals to 1
#     |x, y, z> to |x, z, y>
function cswap(state::QuantumState, ctrl::Int, x::Int, y::Int)
    newstate = Vector{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, ctrl)
            newstate[i] = state[i]
        else
            mapped = swapbit(basis, x, y)
            newstate[i] = state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

# Performs the following mapping:
# |00> to |00>
# |01> to (0.5 + 0.5i) * |01> + (0.5 - 0.5i) * |10>
# |10> to (0.5 - 0.5i) * |01> + (0.5 + 0.5i) * |10>
# |11> to |11>
# This more generally translates to the following, less efficient
# implementation:
# |x, y> to (0.5 + 0.5i) * |x, y> + (0.5 - 0.5i) * |y, x>
function sqrtswap(state::QuantumState, x::Int, y::Int)
    newstate = Vector{Complex{Float64}}(length(state))

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

# Given theta, it performs the following mapping:
# |0> to |0>
# |1> to e^(i * theta) * |1>
function phaseshift(state::QuantumState, bit::Int, theta::Float64)
    newstate = Vector{Complex{Float64}}(length(state))
    bit = state.bits - bit
    mul = exp(im * theta)

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, bit)
            newstate[i] = state[i]
        else
            newstate[i] = mul * state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

function cphase(state::QuantumState, control::Int, bit::Int, theta::Float64)
    newstate = Vector{Complex{Float64}}(length(state))
    control = state.bits - bit
    bit = state.bits - bit
    mul = exp(im * theta)

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, bit) || iszerobit(basis, control)
            newstate[i] = state[i]
        else
            newstate[i] = mul * state[i]
        end
    end

    return QuantumState(newstate, state.bits)
end

paulix(state::QuantumState, bit::Int) = not(state, bit)

# Performs the following mapping:
# |0> to i|1>
# |1> to -i|0>
function pauliy(state::QuantumState, bit::Int)
    newstate = Vector{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, bit)
            mapped = setbit(basis, bit)
            newstate[i] =  im * state[mapped + 1]
        else
            mapped = clearbit(basis, bit)
            newstate[i] = -im * state[mapped + 1]
        end
    end

    return QuantumState(newstate, state.bits)
end

# A special case of the phaseshift gate. Performs the following mapping:
# |0> to |0>
# |1> to -|1>
function pauliz(state::QuantumState, bit::Int)
    newstate = Vector{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, bit)
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

    # First swap the bits to the first N positions
    for i = 1:length(bits)
        newstate = swap(newstate, i, bits[i])

        # Track swaps that would effect later qubits. An example:
        # Given a 3 bit state with an operation defined on bits 3 and 1 in that order,
        # let the bits be the x, y and x, their indices 1 2 and 3. We want to achieve
        # the following operation: x, y, z -> z, x, y
        # Step 0: Initial state: x, y, z [remaining bits: 3, 1]
        # Step 1: Swap the bit at position 1 with the first bit of the operation:
        #         x, y, z -> z, y, x [remaining bits: 1]
        # Step 2: Update the remaing bits. Replace every reference to position 1
        #         with the position it's being swapped to: [remaining bits: 3]
        # Step 3: Proceed to Step 1 with the updated remaining list
        for j = i + 1:length(bits)
            if bits[j] == i
                bits[j] = bits[i]
            end
        end
    end

    # Perform the unitary operation
    if length(state) == size(matrix, 1)
        newstate = QuantumState(matrix * newstate.vector, newstate.bits)
    else
        gatematrix = kron(matrix, speye(2^(state.bits - length(bits))))
        newstate = QuantumState(gatematrix * newstate.vector, newstate.bits)
    end

    # The bits variable stores the actual operations that were performed to achieve
    # the current bit order. Performing these in reverse order restores the correct
    # qubit order
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

# Chooses one of the given gates, based on the decimal representation of the classical bit state given in
# the bits argument. It performs the following mapping:
# |x> to gates[decimal(bits) + 1](|x>)
# Where the +1 adjusts the result to Julia's 1-based indexing system
function choose1(state::QuantumState, bits::Vector{Int}, gates::Vector{Vector{Tuple{Function,Vector{Any}}}})
    index = todecimal(bits) + 1
    chain = gates[index]

    newstate::QuantumState = chain[1][1](state, chain[1][2]...)
    for i = 2:length(chain)
        newstate = chain[i][1](newstate, chain[i][2]...)
    end
    return newstate
end

# Selects a random basis vector based on the probabilities given in the argument.
# It uses the same algorithm as the one used in Fitness Proportionate Selection.
# TODO(gustorn): See if the stochastic acceptance variant performs better for large
#                quantum states
function randbasis(probabilities::Vector{T}) where {T <: Real}
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

# Measures the quantum state, collapsing it and returning both the resulting
# posterior state (with the appropriate basis vector having a probability of 1.0)
# and the result of the measurement.
function measure(state::QuantumState)
    probabilities = map(abs2, state.vector)::Vector{Float64}
    basis = randbasis(probabilities)

    # Measure the quantum state
    measurement = reverse(digits(basis, 2, state.bits))

    # Collapse the state
    newstate = fill(0.0 + 0im, length(state))
    newstate[basis + 1] = 1.0 + 0im
    posterior = QuantumState(newstate, state.bits)

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

    # First collect the coefficients for the probabilities of the qubit being |0>
    prob0 = 0.0
    for i = 1:length(state)
        basis = i - 1
        if iszerobit(basis, bit)
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
    n = sqrt(sum(abs2, posterior))
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

    map!(b -> state.bits - b, bits, bits)
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
    n = sqrt(sum(abs2, posterior))
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


