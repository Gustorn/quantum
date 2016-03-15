module Gates

using Iterators.zip

using QSpice.State
using QSpice.Util

export identity, hadamard, rotation, not, cnot, swap,
       phase_shift, pauli_x, pauli_y, pauli_z,
       measure, partial_measure, probe

function identity(state::QuantumState)
    return copy(state)
end

function hadamard(state::QuantumState, bit::Int)
    new_state = fill(0.0 + 0im, length(state))
    bit = state.bits - bit
    invsqrt2 = 1.0 / sqrt(2)

    for i = 1:length(state)
        b  = i - 1
        b0 = clear_bit(b, bit)
        b1 = set_bit(b, bit)

        new_state[b0 + 1] += invsqrt2 * state[i]
        if is_zero(b, bit)
            new_state[b1 + 1] += invsqrt2 * state[i]
        else
            new_state[b1 + 1] -= invsqrt2 * state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

function not(state::QuantumState, bit::Int)
    new_state = fill(0.0 + 0im, length(state))

    bit = state.bits - bit
    for i = 1:length(state)
        b = i - 1
        new_b = flip_bit(b, bit)
        new_state[new_b + 1] += state[i]
    end

    return QuantumState(new_state, state.bits)
end

function cnot(state::QuantumState, ctrl::Int, flip::Int)
    new_state = fill(0.0 + 0im, length(state))

    ctrl = state.bits - ctrl
    flip = state.bits - flip

    for i = 1:length(state)
        b = i - 1
        if is_zero(b, ctrl)
            new_state[i] += state[i]
        else
            new_b = flip_bit(b, flip)
            new_state[new_b + 1] += state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

function swap(state::QuantumState, x::Int, y::Int)
    new_state = fill(0.0 + 0im, length(state))

    x = state.bits - x
    y = state.bits - y

    for i = 1:length(state)
        b = i - 1
        b_new = swap_bit(b, x, y)
        new_state[b_new + 1] = state[i]
    end

    return QuantumState(new_state, state.bits)
end

function phase_shift(state::QuantumState, bit::Int, theta::Float64)
    new_state = fill(0.0 + 0im, length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        b = i - 1
        if is_zero(b, bit)
            new_state[i] += state[i]
        else
            new_state[i] += exp(im * theta) * state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

pauli_x(state::QuantumState, bit::Int) = not(state, bit)

function pauli_y(state::QuantumState, bit::Int)
    new_state = fill(0.0 + 0im, length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        b = i - 1
        if is_zero(b, bit)
            b_new = set_bit(b, bit)
            new_state[i] += im * state[b_new + 1]
        else
            b_new = clear_bit(b, bit)
            new_state[i] += -im * state[b_new + 1]
        end
    end

    return QuantumState(new_state, state.bits)
end

function pauli_z(state::QuantumState, bit::Int)
    new_state = fill(0.0 + 0im, length(state))
    bit = state.bits - bit

    for i = 1:length(state)
        b = i - 1
        if is_zero(b, bit)
            new_state[i] += state[i]
        else
            new_state[i] -= state[i]
        end
    end

    return QuantumState(new_state, state.bits)
end

# This acts as a terminator right now. It could return a pure qubit
# state that matches the measurement, but generally this doesn't really
# have any benefits
function measure(state::QuantumState)
    probabilities = map(abs2, state.vector)
    reverse(digits(weighed_random_state(probabilities), 2, state.bits))
end

# Performs partial measurement on the given bit of the quantum state
# in the computational basis
function partial_measure(state::QuantumState, bit::Int)
    new_state = fill(0.0 + 0im, 2^(state.bits - 1))
    prob0 = 0.0
    prob1 = 0.0

    bit = state.bits - bit

    # First collect the coefficients for the probabilities of the
    # |0> and |1> states respectively
    for b = 0:length(state) - 1
        if is_zero(b, bit)
            prob0 += abs2(state[b + 1])
        else
            prob1 += abs2(state[b + 1])
        end
    end

    # And measure the bit
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

    return (measurement, QuantumState(new_state, state.bits - 1))
end

function partial_measure(state::QuantumState, x0::Int, x1::Int, xs::Int...)
    # First construct the array of bits we want to measure.
    bits = [x0, x1]
    append!(bits, [xs...])
    map!(b -> state.bits - b, bits)
    sort!(bits)

    num_measured  = length(bits)
    num_posterior = state.bits - num_measured
    probabilities = fill(0.0, 2^num_measured)
    new_state = fill(0.0 + 0im, 2^num_posterior)

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
        for bit in bits
            cb = get_bit(basis, bit)
            partial_index = partial_index | (cb << shift_by)
            shift_by += 1
        end

        # Add |a|^2 to the probability of the appropriate basis
        probabilities[partial_index + 1] += abs2(state[basis + 1])
    end

    # Measure it, keeping the same format as the regular `measure` function
    measurement = reverse(digits(weighed_random_state(probabilities), 2, num_measured))


    # The measurement returns the result in 1:bit[0], length:bit[length] order, while
    # the original bits array stores them in 1:bit[length] length: bit[0]
    bits_rev = reverse(bits)

    # Now build up the posterior state
    for basis in 0:length(state) - 1
        # Check if the current basis vector matches the measurement at the specified bits
        state_matches = true
        for (i, bit) in zip(bits_rev, measurement)
            if get_bit(basis, i) != bit
                state_matches = false
            end
        end

        # We only need to bother with the coefficients if the measurement matches the
        # partial basis vector
        if state_matches
            # Now we construct the basis vector of the new state by removing the part we
            # measured. It's important that this operation is performed in LSB->MSB order:
            # this way we can consistently lower the index of the bit we want to remove,
            # because the next bit is guaranteed to be higher than previous one
            removed_bits = 0
            new_index = basis
            for bit in bits
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

    return (measurement, QuantumState(new_state, state.bits - length(bits)))
end

# Prints the current state of the qubit in terms of basis vectors, coefficients
# and probabilities. Useful for debugging and visualization. Guaranteed not to disturb
# the quantum state
function probe(state::QuantumState, name::AbstractString)
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


