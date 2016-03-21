module Gates

using Iterators

using QSpice.State
using QSpice.Util

export identity, hadamard, not, cnot, swap, sqrt_swap,
       phase_shift, pauli_x, pauli_y, pauli_z,
       ccnot, cswap, unitary,
       measure, partial_measure, probe

function identity(state::QuantumState)
    return copy(state)
end

function hadamard(state::QuantumState, bit::Int)
    new_state = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit
    invsqrt2 = 1.0 / sqrt(2)

    for i = 1:length(state)
        basis  = i - 1
        mapped0 = clear_bit(basis, bit)
        mapped1 = set_bit(basis, bit)

        new_state[i] = invsqrt2 * state[mapped0 + 1]
        if is_zero(basis, bit)
            new_state[i] += invsqrt2 * state[mapped1 + 1]
        else
            new_state[i] -= invsqrt2 * state[mapped1 + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function not(state::QuantumState, bit::Int)
    new_state = Array{Complex{Float64}}(length(state))

    bit = state.bits - bit
    for i = 1:length(state)
        basis = i - 1
        mapped = flip_bit(basis, bit)
        new_state[i] = state[mapped + 1]
    end

    return QuantumState(new_state, state.bits)
end

function cnot(state::QuantumState, ctrl::Int, flip::Int)
    new_state = Array{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    flip = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if is_zero(basis, ctrl)
            new_state[i] = state[i]
        else
            mapped = flip_bit(basis, flip)
            new_state[i] = state[mapped + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function ccnot(state::QuantumState, ctrl1::Int, ctrl2::Int, flip::Int)
    new_state = Array{Complex{Float64}}(length(state))

    ctrl1 = state.bits - ctrl1
    ctrl2 = state.bits - ctrl2
    flip  = state.bits - flip

    for i = 1:length(state)
        basis = i - 1
        if !is_zero(basis, ctrl1) && !is_zero(basis, ctrl2)
            mapped = flip_bit(basis, flip)
            new_state[i] = state[mapped + 1]
        else
            new_state[i] = state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

function swap(state::QuantumState, x::Int, y::Int)
    new_state = Array{Complex{Float64}}(length(state))

    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        mapped = swap_bit(basis, x, y)
        new_state[i] = state[mapped + 1]
    end

    return QuantumState(new_state, state.bits)
end

function cswap(state::QuantumState, ctrl::Int, x::Int, y::Int)
    new_state = Array{Complex{Float64}}(length(state))

    ctrl = state.bits - ctrl
    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        basis = i - 1
        if is_zero(basis, ctrl)
            new_state[i] = state[i]
        else
            mapped = swap_bit(basis, x, y)
            new_state[i] = state[mapped + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function sqrt_swap(state::QuantumState, x::Int, y::Int)
    new_state = Array{Complex{Float64}}(length(state))

    x = state.bits - x
    y = state.bits - y

    mul_orig   = 0.5 + 0.5im
    mul_mapped = 0.5 - 0.5im

    for i = 1:length(state)
        basis = i - 1
        mapped = swap_bit(basis, x, y)
        if basis == mapped
            new_state[i] = state[i]
        else
            new_state[i] = mul_orig * state[i] + mul_mapped * state[mapped + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function phase_shift(state::QuantumState, bit::Int, theta::Float64)
    new_state = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit
    etheta = exp(im * theta)

    for i = 1:length(state)
        basis = i - 1
        if is_zero(basis, bit)
            new_state[i] = state[i]
        else
            new_state[i] = etheta * state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

pauli_x(state::QuantumState, bit::Int) = not(state, bit)

function pauli_y(state::QuantumState, bit::Int)
    new_state = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if is_zero(basis, bit)
            mapped = set_bit(basis, bit)
            new_state[i] =  im * state[mapped + 1]
        else
            mapped = clear_bit(basis, bit)
            new_state[i] = -im * state[mapped + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function pauli_z(state::QuantumState, bit::Int)
    new_state = Array{Complex{Float64}}(length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        basis = i - 1
        if is_zero(basis, bit)
            new_state[i] =  state[i]
        else
            new_state[i] = -state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

function unitary(state::QuantumState, matrix::Array)
    if length(state.vector) != size(matrix, 1) || length(state.vector) != size(matrix, 2)
        error("The unitary matrix dimensions don't match the quantum state vector's");
    end

    return QuantumState(matrix * state.vector, state.bits)
end

function unitary(state::QuantumState, matrix::Array, bit0::Int, bits::Int...)
    bits = [bit0, bits...]
    map!(b -> state.bits - b, bits)
    mapped_state_size = 2^length(bits)

    # In case the matrix is a perfect match, a straight-up multiply is much faster
    if length(state) == size(matrix, 1)
        return unitary(state, matrix)
    end

    if mapped_state_size != size(matrix, 1) || mapped_state_size != size(matrix, 2)
        error("The unitary matrix does not map enough bits")
    end

    # Cannot take the kronecker products of a single qubit state
    if length(bits) == 1
        mapped_bases = Vector[QUBIT_0.vector, QUBIT_1.vector]
    else
        # Generate the basis vectors in the order of: |00..00>, |00..01>, |00..10>, etc...
        mapped_bases = [kron([get_bit(basis, qubit) == 0 ? QUBIT_0.vector : QUBIT_1.vector
                              for qubit in length(bits)-1:-1:0]...)
                        for basis in 0:mapped_state_size - 1]
    end
    map!(x -> matrix * x, mapped_bases)

    new_state = fill(0.0 + 0im, length(state))
    for i = 1:length(state)
        basis = i - 1

        # First we construct the index for the mapped basis vectors
        partial_index = 0
        shift_by = length(bits) - 1
        for bit in bits
            cb = get_bit(basis, bit)
            partial_index = partial_index | (cb << shift_by)
            shift_by -= 1
        end
        partial_index += 1

        # Then for each mapped basis vector we collect the appropriate coefficients
        # from the original states
        for j=1:length(mapped_bases[partial_index]), k = 1:length(state)
            mapped = j - 1
            try_match = k - 1

            # Check if the current (original) basis matches the mapped basis vector
            # at the appropriate qubits
            is_matching = true
            for (mi, bit) in enumerate(bits)
                if get_bit(try_match, bit) != get_bit(mapped, length(bits) - mi)
                    is_matching = false;
                    break
                end
            end

            if is_matching
                # If it matched then collect the coefficients, if the rest of the basis vectors match
                inter_basis = reduce((a, x) -> clear_bit(a, x), basis, bits)
                inter_match = reduce((a, x) -> clear_bit(a, x), try_match, bits)

                if inter_basis == inter_match
                    new_state[i] += mapped_bases[partial_index][mapped + 1] * state[try_match + 1]
                end
            end
        end
    end

    return QuantumState(new_state, state.bits)
end

# This acts as a terminator right now. It could return a pure qubit
# state that matches the measurement, but generally this doesn't really
# have any benefits
function measure(state::QuantumState)
    probabilities = map(abs2, state.vector)
    return reverse(digits(weighed_random_state(probabilities), 2, state.bits))
end

# Performs partial measurement on the given bit of the quantum state
# in the computational basis
function partial_measure(state::QuantumState, bit::Int)
    if state.bits == 1
        return (measure(state), EMPTY_STATE)
    end

    new_state::Vector{Complex{Float64}} = fill(0.0 + 0im, 2^(state.bits - 1))

    # We can assume the qubit is normalized, so we only need to figure
    # out the probability of getting a |0>. If that didn't happen it's
    # guaranteed to be |1>
    prob0 = 0.0

    bit = state.bits - bit

    # First collect the coefficients for the probabilities of the
    # |0> and |1> states respectively
    for b = 0:length(state) - 1
        if is_zero(b, bit)
            prob0 += abs2(state[b + 1])
        end
    end

    # And measure the bit.
    measurement = rand() <= prob0 ? 0 : 1

    # Next collect the coefficients for the basis vectors of the
    # posterior states
    for b = 0:length(state) - 1
        if get_bit(b, bit) == measurement
            new_basis = shift_range_down(b, bit)
            new_state[new_basis + 1] += state[b + 1]
        end
    end

    # Now normalize the new states so that the sum of coefficients is equal
    # to 1
    # TODO(gustorn): see if the error from FP precision has any visible
    # effect on the posterior states
    n = sqrt(reduce((accum, x) -> accum + abs2(x), 0.0, new_state))
    map!(x -> x / n, new_state)

    return ([measurement], QuantumState(new_state, state.bits - 1))
end

function partial_measure(state::QuantumState, x0::Int, x1::Int, xs::Int...)
    # First construct the array of bits we want to measure.
    measured_bits = [x0, x1, xs...]
    measured_bits = unique(measured_bits)

    if length(measured_bits) == state.bits
        return (measure(state), EMPTY_STATE)
    end

    # Next perform the usual mapping: bit0 is physically represented by a higher
    # bit in the index than bitN. Example:
    # Given the state: |1000>, the index representing qubit0, which in case is |1>,
    # is 0b1000, hence the need for the subtraction
    map!(b -> state.bits - b, measured_bits)

    # And finally sort the bits, so the lowest (phyiscal) / highest (real) bits
    # come first.
    sort!(measured_bits)

    num_measured  = length(measured_bits)
    num_posterior = state.bits - num_measured
    probabilities::Vector{Float64} = fill(0.0, 2^num_measured)
    new_state::Vector{Complex{Float64}} = fill(0.0 + 0im, 2^num_posterior)

    # First build up the probabilities for the bits we're going to measure
    for basis in 0:length(state) - 1
        # First construct the index of the basis vector in the
        # new, partial state. This is done by incrementally copying
        # the specific bits of the old basis vector into a new one
        # Example: consider the basis vector |100> where want to copy
        # bits 1 and 3 (|1> and |0> respectively). First we copy |0> into
        # the first position of the new basis vector, giving us |x0>, where x
        # is not yet filled in. Next we copy the |1>, yielding the desired
        # |10> state
        partial_index = 0
        shift_by = 0
        for bit in measured_bits
            cb = get_bit(basis, bit)
            partial_index = partial_index | (cb << shift_by)
            shift_by += 1
        end

        # Add |a|^2 to the probability of the appropriate basis
        probabilities[partial_index + 1] += abs2(state[basis + 1])
    end

    # Measure it. Important note: this is the reverse of the order returned by
    # measure, but it makes further calculations easier. Don't forget to reverse
    # it before returning it.
    measurement = digits(weighed_random_state(probabilities), 2, num_measured)

    # Now build up the posterior state
    for basis in 0:length(state) - 1
        # Check if the current basis vector matches the measurement at the specified bits
        state_matches = true
        @itr for (i, bit) in zip(measured_bits, measurement)
            if get_bit(basis, i) != bit
                state_matches = false
                break
            end
        end

        # We only need to bother with the coefficients if the measurement matches the
        # partial basis vector
        if state_matches
            # Now we construct the basis vector of the new state by removing the part we
            # measured. It's important that this operation is performed in LSB->MSB order:
            # this way we can consistently lower the index of the bit we want to remove,
            # because the next bit is guaranteed to be higher than previous one. Example:
            # Given the state |1010>, if we measured qubits 2 and 3 (|01> in this example),
            # then the algorithm performs the following:
            # 1. First remove the first (lowest physical/highest actual) bit at index 2,
            #    highlighted: |10(1)0>
            # 2. This yields the following state: |100>. The next index we want to remove is 3,
            #    but note how that index shifted down by one place. With a more visual representation,
            #    step 1 performed the following: |abcd> -> |abd>, or with indices: |4321> -> |431>
            # 3. Decrement the index we want to remove at, and proceed with step 1
            removed_bits = 0
            new_index = basis
            for bit in measured_bits
                new_index = shift_range_down(new_index, bit - removed_bits)
                removed_bits += 1
            end

            # Add the coefficients from the old state to the new one
            new_state[new_index + 1] += state[basis + 1]
        end
    end

    # Normalize the new state so that sum(|coeff|^2) == 1.0
    # The same concerns about FP precision apply here just as
    # they do in the single-bit version of partial measurement
    n = sqrt(reduce((accum, x) -> accum + abs2(x), 0.0, new_state))
    map!(x -> x / n, new_state)

    # Match the bit order of the regular `measure` function
    reverse!(measurement)
    return (measurement, QuantumState(new_state, state.bits - length(measured_bits)))
end

# Prints the current state of the qubit in terms of basis vectors, coefficients
# and probabilities. Useful for debugging and visualization. Guaranteed not to disturb
# the quantum state
function probe{T<:AbstractString}(state::QuantumState, name::T)
    println(name)
    print_bases(copy(state))
    println()
end

# The same as probe, but without the need for an explicit name
function probe(state::QuantumState)
    print_bases(copy(state))
    println()
end

end


