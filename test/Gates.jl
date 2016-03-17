include("Common.jl")

using FactCheck
using QSpice.Gates, QSpice.State

const RUNS_PER_STATE = 100
const QUBITS = [QUBIT_0, QUBIT_1]

random_state(nbits::Int) = from_states([rand(QUBITS) for _ in 1:nbits]...)

function random_test_range(fn, bits_from, bits_to)
    for i = bits_from:bits_to
        for _=1:RUNS_PER_STATE
            test_state = random_state(i)
            fn(test_state)
        end
    end
end

facts("Quantum state creation sanity tests") do
    context("Single qubit state") do
        @fact QUBIT_0 --> roughly([1.0 + 0im, 0.0 + 0im])
        @fact QUBIT_1 --> roughly([0.0 + 0im, 1.0 + 0im])
    end

    context("Basis vectors for N qubit states") do
        random_test_range(2, 10) do test_state
            @fact countnz(test_state.vector) --> 1
        end
    end
end

facts("Hadamard gate") do
    context("Single qubit mapping") do
        @fact hadamard(QUBIT_0, 1) --> roughly((1.0 / sqrt(2)) .* [1.0 + 0im,  1])
        @fact hadamard(QUBIT_1, 1) --> roughly((1.0 / sqrt(2)) .* [1.0 + 0im, -1])
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        isq = 1.0 / sqrt(2)
        @fact hadamard(s00, 1) --> roughly(isq .* [1.0 + 0im,  0,  1,  0])
        @fact hadamard(s00, 2) --> roughly(isq .* [1.0 + 0im,  1,  0,  0])

        @fact hadamard(s01, 1) --> roughly(isq .* [0.0 + 0im,  1,  0,  1])
        @fact hadamard(s01, 2) --> roughly(isq .* [1.0 + 0im, -1,  0,  0])

        @fact hadamard(s10, 1) --> roughly(isq .* [1.0 + 0im,  0, -1,  0])
        @fact hadamard(s10, 2) --> roughly(isq .* [0.0 + 0im,  0,  1,  1])

        @fact hadamard(s11, 1) --> roughly(isq .* [0.0 + 0im,  1,  0, -1])
        @fact hadamard(s11, 2) --> roughly(isq .* [0.0 + 0im,  0,  1, -1])
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            @fact hadamard(hadamard(test_state, x), x) --> roughly(test_state)
        end
    end
end

facts("Not gate") do
    context("Single qubit mapping") do
        @fact Gates.not(QUBIT_0, 1) --> roughly(QUBIT_1)
        @fact Gates.not(QUBIT_1, 1) --> roughly(QUBIT_0)
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact Gates.not(s00, 1) --> roughly(s10)
        @fact Gates.not(s00, 2) --> roughly(s01)

        @fact Gates.not(s01, 1) --> roughly(s11)
        @fact Gates.not(s01, 2) --> roughly(s00)

        @fact Gates.not(s10, 1) --> roughly(s00)
        @fact Gates.not(s10, 2) --> roughly(s11)

        @fact Gates.not(s11, 1) --> roughly(s01)
        @fact Gates.not(s11, 2) --> roughly(s10)
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            @fact Gates.not(Gates.not(test_state, x), x) --> roughly(test_state)
        end
    end
end


facts("CNOT gate") do
    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact cnot(s00, 1, 2) --> roughly(s00)
        @fact cnot(s00, 2, 1) --> roughly(s00)

        @fact cnot(s01, 1, 2) --> roughly(s01)
        @fact cnot(s01, 2, 1) --> roughly(s11)

        @fact cnot(s10, 1, 2) --> roughly(s11)
        @fact cnot(s10, 2, 1) --> roughly(s10)

        @fact cnot(s11, 1, 2) --> roughly(s10)
        @fact cnot(s11, 2, 1) --> roughly(s01)
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits - 1)
            y = rand(x + 1:test_state.bits)
            @fact cnot(cnot(test_state, x, y), x, y) --> roughly(test_state)
        end
    end
end

facts("CCNOT gate") do
    context("Three qubit mapping") do
        s000 = from_states(QUBIT_0, QUBIT_0, QUBIT_0)
        s001 = from_states(QUBIT_0, QUBIT_0, QUBIT_1)
        s010 = from_states(QUBIT_0, QUBIT_1, QUBIT_0)
        s011 = from_states(QUBIT_0, QUBIT_1, QUBIT_1)

        s100 = from_states(QUBIT_1, QUBIT_0, QUBIT_0)
        s101 = from_states(QUBIT_1, QUBIT_0, QUBIT_1)
        s110 = from_states(QUBIT_1, QUBIT_1, QUBIT_0)
        s111 = from_states(QUBIT_1, QUBIT_1, QUBIT_1)

        @fact ccnot(s000, 1, 2, 3) --> roughly(s000)
        @fact ccnot(s001, 1, 2, 3) --> roughly(s001)
        @fact ccnot(s010, 1, 2, 3) --> roughly(s010)
        @fact ccnot(s011, 1, 2, 3) --> roughly(s011)
        @fact ccnot(s100, 1, 2, 3) --> roughly(s100)
        @fact ccnot(s101, 1, 2, 3) --> roughly(s101)
        @fact ccnot(s110, 1, 2, 3) --> roughly(s111)
        @fact ccnot(s111, 1, 2, 3) --> roughly(s110)

        @fact ccnot(s011, 2, 3, 1) --> roughly(s111)
        @fact ccnot(s101, 1, 3, 2) --> roughly(s111)
    end

    context("Unitarity on N qubit states") do
        random_test_range(3, 10) do test_state
            x = rand(1:test_state.bits - 2)
            y = rand(x + 1:test_state.bits - 1)
            z = rand(y + 1:test_state.bits)
            @fact ccnot(ccnot(test_state, x, y, z), x, y, z) --> roughly(test_state)
        end
    end
end


facts("Swap gate") do
    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact swap(s00, 1, 2) --> roughly(s00)
        @fact swap(s00, 2, 1) --> roughly(s00)

        @fact swap(s01, 1, 2) --> roughly(s10)
        @fact swap(s01, 2, 1) --> roughly(s10)

        @fact swap(s10, 1, 2) --> roughly(s01)
        @fact swap(s10, 2, 1) --> roughly(s01)

        @fact swap(s11, 1, 2) --> roughly(s11)
        @fact swap(s11, 2, 1) --> roughly(s11)
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            y = rand(1:test_state.bits)
            @fact swap(swap(test_state, x, y), y, x) --> roughly(test_state)
        end
    end
end

facts("CSWAP gate") do
    context("Three qubit mapping") do
        s000 = from_states(QUBIT_0, QUBIT_0, QUBIT_0)
        s001 = from_states(QUBIT_0, QUBIT_0, QUBIT_1)
        s010 = from_states(QUBIT_0, QUBIT_1, QUBIT_0)
        s011 = from_states(QUBIT_0, QUBIT_1, QUBIT_1)

        s100 = from_states(QUBIT_1, QUBIT_0, QUBIT_0)
        s101 = from_states(QUBIT_1, QUBIT_0, QUBIT_1)
        s110 = from_states(QUBIT_1, QUBIT_1, QUBIT_0)
        s111 = from_states(QUBIT_1, QUBIT_1, QUBIT_1)

        @fact cswap(s000, 1, 2, 3) --> roughly(s000)
        @fact cswap(s001, 1, 2, 3) --> roughly(s001)
        @fact cswap(s010, 1, 2, 3) --> roughly(s010)
        @fact cswap(s011, 1, 2, 3) --> roughly(s011)
        @fact cswap(s100, 1, 2, 3) --> roughly(s100)
        @fact cswap(s101, 1, 2, 3) --> roughly(s110)
        @fact cswap(s110, 1, 2, 3) --> roughly(s101)
        @fact cswap(s111, 1, 2, 3) --> roughly(s111)

        @fact cswap(s011, 2, 1, 3) --> roughly(s110)
        @fact cswap(s101, 3, 1, 2) --> roughly(s011)
    end

    context("Unitarity on N qubit states") do
        random_test_range(3, 10) do test_state
            x = rand(1:test_state.bits - 2)
            y = rand(x + 1:test_state.bits - 1)
            z = rand(y + 1:test_state.bits)
            @fact cswap(cswap(test_state, x, y, z), x, z, y) --> roughly(test_state)
        end
    end
end

facts("Square root swap gate") do
    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact sqrt_swap(s00, 1, 2) --> roughly(s00)
        @fact sqrt_swap(s00, 2, 1) --> roughly(s00)

        @fact sqrt_swap(s01, 1, 2) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0])
        @fact sqrt_swap(s01, 2, 1) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0])

        @fact sqrt_swap(s10, 1, 2) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0])
        @fact sqrt_swap(s10, 2, 1) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0])

        @fact sqrt_swap(s11, 1, 2) --> roughly(s11)
        @fact sqrt_swap(s11, 2, 1) --> roughly(s11)
    end

    context("Three qubit mapping") do
        s000 = from_states(QUBIT_0, QUBIT_0, QUBIT_0)
        s001 = from_states(QUBIT_0, QUBIT_0, QUBIT_1)
        s010 = from_states(QUBIT_0, QUBIT_1, QUBIT_0)
        s011 = from_states(QUBIT_0, QUBIT_1, QUBIT_1)

        s100 = from_states(QUBIT_1, QUBIT_0, QUBIT_0)
        s101 = from_states(QUBIT_1, QUBIT_0, QUBIT_1)
        s110 = from_states(QUBIT_1, QUBIT_1, QUBIT_0)
        s111 = from_states(QUBIT_1, QUBIT_1, QUBIT_1)
        @fact sqrt_swap(s000, 1, 2) --> roughly(s000)
        @fact sqrt_swap(s000, 2, 1) --> roughly(s000)
        @fact sqrt_swap(s000, 1, 3) --> roughly(s000)
        @fact sqrt_swap(s000, 3, 1) --> roughly(s000)
        @fact sqrt_swap(s000, 2, 3) --> roughly(s000)
        @fact sqrt_swap(s000, 3, 2) --> roughly(s000)

        @fact sqrt_swap(s001, 1, 2) --> roughly(s001)
        @fact sqrt_swap(s001, 2, 1) --> roughly(s001)
        @fact sqrt_swap(s001, 1, 3) --> roughly([0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrt_swap(s001, 3, 1) --> roughly([0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrt_swap(s001, 2, 3) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0, 0, 0, 0, 0])
        @fact sqrt_swap(s001, 3, 2) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0, 0, 0, 0, 0])

        @fact sqrt_swap(s010, 1, 2) --> roughly([0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrt_swap(s010, 2, 1) --> roughly([0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrt_swap(s010, 1, 3) --> roughly(s010)
        @fact sqrt_swap(s010, 3, 1) --> roughly(s010)
        @fact sqrt_swap(s010, 2, 3) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0, 0, 0, 0, 0])
        @fact sqrt_swap(s010, 3, 2) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0, 0, 0, 0, 0])

        @fact sqrt_swap(s011, 1, 2) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0])
        @fact sqrt_swap(s011, 2, 1) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0])
        @fact sqrt_swap(s011, 1, 3) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0])
        @fact sqrt_swap(s011, 3, 1) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0])
        @fact sqrt_swap(s011, 2, 3) --> roughly(s011)
        @fact sqrt_swap(s011, 3, 2) --> roughly(s011)

        @fact sqrt_swap(s100, 1, 2) --> roughly([0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrt_swap(s100, 2, 1) --> roughly([0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrt_swap(s100, 1, 3) --> roughly([0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrt_swap(s100, 3, 1) --> roughly([0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrt_swap(s100, 2, 3) --> roughly(s100)
        @fact sqrt_swap(s100, 3, 2) --> roughly(s100)

        @fact sqrt_swap(s101, 1, 2) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0])
        @fact sqrt_swap(s101, 2, 1) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0])
        @fact sqrt_swap(s101, 1, 3) --> roughly(s101)
        @fact sqrt_swap(s101, 3, 1) --> roughly(s101)
        @fact sqrt_swap(s101, 2, 3) --> roughly([0, 0, 0, 0, 0, 0.5 + 0.5im, 0.5 - 0.5im, 0])
        @fact sqrt_swap(s101, 3, 2) --> roughly([0, 0, 0, 0, 0, 0.5 + 0.5im, 0.5 - 0.5im, 0])

        @fact sqrt_swap(s110, 1, 2) --> roughly(s110)
        @fact sqrt_swap(s110, 2, 1) --> roughly(s110)
        @fact sqrt_swap(s110, 1, 3) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0])
        @fact sqrt_swap(s110, 3, 1) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0])
        @fact sqrt_swap(s110, 2, 3) --> roughly([0, 0, 0, 0, 0, 0.5 - 0.5im, 0.5 + 0.5im, 0])
        @fact sqrt_swap(s110, 3, 2) --> roughly([0, 0, 0, 0, 0, 0.5 - 0.5im, 0.5 + 0.5im, 0])

        @fact sqrt_swap(s111, 1, 2) --> roughly(s111)
        @fact sqrt_swap(s111, 2, 1) --> roughly(s111)
        @fact sqrt_swap(s111, 1, 3) --> roughly(s111)
        @fact sqrt_swap(s111, 3, 1) --> roughly(s111)
        @fact sqrt_swap(s111, 2, 3) --> roughly(s111)
        @fact sqrt_swap(s111, 3, 2) --> roughly(s111)
    end
end

facts("Phase shift gate") do
    context("Single qubit mapping") do
        @fact phase_shift(QUBIT_0, 1, pi / 8) --> roughly(QUBIT_0)
        @fact phase_shift(QUBIT_1, 1, pi / 8) --> roughly([0.0 + 0im, exp(im * pi / 8)])
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)

        @fact phase_shift(s00, 1, pi / 8) --> roughly(s00)
        @fact phase_shift(s00, 2, pi / 8) --> roughly(s00)

        @fact phase_shift(s01, 1, pi / 8) --> roughly(s01)
        @fact phase_shift(s01, 2, pi / 8) --> roughly([0, exp(im * pi / 8), 0, 0])

        @fact phase_shift(s10, 1, pi / 8) --> roughly([0, 0, exp(im * pi / 8), 0])
        @fact phase_shift(s10, 2, pi / 8) --> roughly(s10)

        @fact phase_shift(s11, 1, pi / 8) --> roughly([0, 0, 0, exp(im * pi / 8)])
        @fact phase_shift(s11, 2, pi / 8) --> roughly([0, 0, 0, exp(im * pi / 8)])
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            theta = rand(1.0:pi)
            @fact phase_shift(phase_shift(test_state, x, theta), x, -theta) --> roughly(test_state)
        end
    end
end

facts("Pauli-X gate") do
    context("Single qubit mapping") do
        @fact pauli_x(QUBIT_0, 1) --> roughly(QUBIT_1)
        @fact pauli_x(QUBIT_1, 1) --> roughly(QUBIT_0)
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact pauli_x(s00, 1) --> roughly(s10)
        @fact pauli_x(s00, 2) --> roughly(s01)

        @fact pauli_x(s01, 1) --> roughly(s11)
        @fact pauli_x(s01, 2) --> roughly(s00)

        @fact pauli_x(s10, 1) --> roughly(s00)
        @fact pauli_x(s10, 2) --> roughly(s11)

        @fact pauli_x(s11, 1) --> roughly(s01)
        @fact pauli_x(s11, 2) --> roughly(s10)
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            @fact pauli_x(pauli_x(test_state, x), x) --> roughly(test_state)
        end
    end
end

facts("Pauli-Y gate") do
    context("Single qubit mapping") do
        @fact pauli_y(QUBIT_0, 1) --> roughly([0.0, -im])
        @fact pauli_y(QUBIT_1, 1) --> roughly([im, 0.0])
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact pauli_y(s00, 1) --> roughly([0.0, 0, -im, 0])
        @fact pauli_y(s00, 2) --> roughly([0.0, -im, 0, 0])

        @fact pauli_y(s01, 1) --> roughly([0.0, 0, 0, -im])
        @fact pauli_y(s01, 2) --> roughly([ im, 0, 0, 0.0])

        @fact pauli_y(s10, 1) --> roughly([ im, 0, 0, 0.0])
        @fact pauli_y(s10, 2) --> roughly([0.0, 0, 0, -im])

        @fact pauli_y(s11, 1) --> roughly([0.0, im, 0, 0])
        @fact pauli_y(s11, 2) --> roughly([0.0, 0, im, 0])
    end
end

facts("Pauli-Z gate") do
    context("Single qubit mapping") do
        @fact pauli_z(QUBIT_0, 1) --> roughly(QUBIT_0)
        @fact pauli_z(QUBIT_1, 1) --> roughly([0.0 + 0im, -1])
    end

    context("Two qubit mapping") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)
        @fact pauli_z(s00, 1) --> roughly(s00)
        @fact pauli_z(s00, 2) --> roughly(s00)

        @fact pauli_z(s01, 1) --> roughly(s01)
        @fact pauli_z(s01, 2) --> roughly([0.0 + 0im, -1, 0, 0])

        @fact pauli_z(s10, 1) --> roughly([0.0 + 0im, 0, -1, 0])
        @fact pauli_z(s10, 2) --> roughly(s10)

        @fact pauli_z(s11, 1) --> roughly([0.0 + 0im, 0, 0, -1])
        @fact pauli_z(s11, 2) --> roughly([0.0 + 0im, 0, 0, -1])
    end

    context("Unitarity on N qubit states") do
        random_test_range(2, 10) do test_state
            x = rand(1:test_state.bits)
            @fact pauli_z(pauli_z(test_state, x), x) --> roughly(test_state)
        end
    end
end

function digits_to_num(darr::Vector{Int})
    result = 0
    offset = length(darr) - 1
    for d in darr
        result |= (d << offset)
        offset -= 1
    end
    return result
end

facts("Measurement") do
    context("Single qubit measurement") do
        @fact measure(QUBIT_0) --> [0]
        @fact measure(QUBIT_1) --> [1]
    end

    context("Equal probability statistical test") do
        s0 = from_vector([0.5, -0.5, 0 - 0.5im, 0 + 0.5im])

        num_runs = 500000
        stats = [0, 0, 0, 0]
        for i=1:num_runs
            n = digits_to_num(measure(s0))
            stats[n + 1] += 1
        end
        @fact stats[1] / num_runs --> roughly(0.25, atol=0.1)
        @fact stats[2] / num_runs --> roughly(0.25, atol=0.1)
        @fact stats[3] / num_runs --> roughly(0.25, atol=0.1)
        @fact stats[4] / num_runs --> roughly(0.25, atol=0.1)
    end

    context("Null probability statistical test") do
        rsqrt8 = 1.0 / sqrt(8)
        s0 = from_vector([rsqrt8 + rsqrt8 * im, rsqrt8 - rsqrt8 * im, 0.0,  1.0 / sqrt(2)])

        num_runs = 500000
        stats = [0, 0, 0, 0]
        for i=1:num_runs
            n = digits_to_num(measure(s0))
            stats[n + 1] += 1
        end
        @fact stats[1] / num_runs --> roughly(0.25, atol=0.1)
        @fact stats[2] / num_runs --> roughly(0.25, atol=0.1)
        @fact stats[3] / num_runs --> 0
        @fact stats[4] / num_runs --> roughly(0.5, atol=0.1)
    end
end

# These tests are not that comprehensive, but they are hard to automate
facts("Partial Measurement") do
    context("Single qubit measurement") do
        @fact partial_measure(QUBIT_0, 1) --> ([0], EMPTY_STATE)
        @fact partial_measure(QUBIT_1, 1) --> ([1], EMPTY_STATE)
    end

    context("Two qubit measurement") do
        s00 = from_states(QUBIT_0, QUBIT_0)
        s01 = from_states(QUBIT_0, QUBIT_1)
        s10 = from_states(QUBIT_1, QUBIT_0)
        s11 = from_states(QUBIT_1, QUBIT_1)

        # State |00>
        m, p = partial_measure(s00, 1)
        @fact m --> [0]
        @fact p --> roughly(QUBIT_0)

        m, p = partial_measure(s00, 2)
        @fact m --> [0]
        @fact p --> roughly(QUBIT_0)

        m, p = partial_measure(s00, 1, 2)
        @fact m --> [0, 0]
        @fact p --> EMPTY_STATE

        # State |01>
        m, p = partial_measure(s01, 1)
        @fact m --> [0]
        @fact p --> roughly(QUBIT_1)

        m, p = partial_measure(s01, 2)
        @fact m --> [1]
        @fact p --> roughly(QUBIT_0)

        m, p = partial_measure(s01, 1, 2)
        @fact m --> [0, 1]
        @fact p --> EMPTY_STATE

        # State |10>
        m, p = partial_measure(s10, 1)
        @fact m --> [1]
        @fact p --> roughly(QUBIT_0)

        m, p = partial_measure(s10, 2)
        @fact m --> [0]
        @fact p --> roughly(QUBIT_1)

        m, p = partial_measure(s10, 1, 2)
        @fact m --> [1, 0]
        @fact p --> EMPTY_STATE

        # State |11>
        m, p = partial_measure(s11, 1)
        @fact m --> [1]
        @fact p --> roughly(QUBIT_1)

        m, p = partial_measure(s11, 2)
        @fact m --> [1]
        @fact p --> roughly(QUBIT_1)

        m, p = partial_measure(s11, 1, 2)
        @fact m --> [1, 1]
        @fact p --> EMPTY_STATE

        # Try states with different probabilities
        h = hadamard(s00, 1)
        for i=1:100
            m, p = partial_measure(h, 1)
            @fact p --> roughly(QUBIT_0)
            m, p = partial_measure(h, 2)
            @fact m --> [0]
            @fact p --> roughly(1.0 / sqrt(2) .* [1.0 + 0im, 1])
        end
    end
end
