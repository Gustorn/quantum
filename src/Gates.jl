module Gates

using QSpice.State
using QSpice.Util

export identity, hadamard, rotation, not, cnot, swap,
       phase_shift, pauli_x, pauli_y, pauli_z,
       measurement, probe

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

# TODO(gustorn): accept 2 different quantum states for the operands
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

# TODO(gustorn): accept 2 different quantum states for the operands
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

# This should somehow interact with entangled states
function measurement(state::QuantumState)
    # TODO(gustorn): proper measurement
    return state
end

function probe(state::QuantumState, name::AbstractString)
    println(name)
    print_bases(copy(state))
    println()
end

end


