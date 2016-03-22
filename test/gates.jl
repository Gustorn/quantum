include("common.jl")

using FactCheck
using QSpice.Gates, QSpice.State

const RUNS_PER_STATE = 100
const QUBITS = [QUBIT0, QUBIT1]

function randomtest(fn, minbits, maxbits)
    for i = minbits:maxbits
        for _=1:RUNS_PER_STATE
            teststate = fromstates([rand(QUBITS) for _ in 1:i]...)
            fn(teststate)
        end
    end
end

facts("Quantum state creation sanity tests") do
    context("Single qubit state") do
        @fact QUBIT0 --> roughly([1.0 + 0im, 0.0 + 0im])
        @fact QUBIT1 --> roughly([0.0 + 0im, 1.0 + 0im])
    end

    context("Basis vectors for N qubit states") do
        randomtest(2, 10) do teststate
            @fact countnz(teststate.vector) --> 1
        end
    end
end

facts("Hadamard gate") do
    context("Single qubit mapping") do
        @fact hadamard(QUBIT0, 1) --> roughly((1.0 / sqrt(2)) .* [1.0 + 0im,  1])
        @fact hadamard(QUBIT1, 1) --> roughly((1.0 / sqrt(2)) .* [1.0 + 0im, -1])
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
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
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            @fact hadamard(hadamard(teststate, x), x) --> roughly(teststate)
        end
    end
end

facts("Not gate") do
    context("Single qubit mapping") do
        @fact Gates.not(QUBIT0, 1) --> roughly(QUBIT1)
        @fact Gates.not(QUBIT1, 1) --> roughly(QUBIT0)
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
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
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            @fact Gates.not(Gates.not(teststate, x), x) --> roughly(teststate)
        end
    end
end


facts("CNOT gate") do
    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
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
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits - 1)
            y = rand(x + 1:teststate.bits)
            @fact cnot(cnot(teststate, x, y), x, y) --> roughly(teststate)
        end
    end
end

facts("CCNOT gate") do
    context("Three qubit mapping") do
        s000 = fromstates(QUBIT0, QUBIT0, QUBIT0)
        s001 = fromstates(QUBIT0, QUBIT0, QUBIT1)
        s010 = fromstates(QUBIT0, QUBIT1, QUBIT0)
        s011 = fromstates(QUBIT0, QUBIT1, QUBIT1)

        s100 = fromstates(QUBIT1, QUBIT0, QUBIT0)
        s101 = fromstates(QUBIT1, QUBIT0, QUBIT1)
        s110 = fromstates(QUBIT1, QUBIT1, QUBIT0)
        s111 = fromstates(QUBIT1, QUBIT1, QUBIT1)

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
        randomtest(3, 10) do teststate
            x = rand(1:teststate.bits - 2)
            y = rand(x + 1:teststate.bits - 1)
            z = rand(y + 1:teststate.bits)
            @fact ccnot(ccnot(teststate, x, y, z), x, y, z) --> roughly(teststate)
        end
    end
end


facts("Swap gate") do
    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
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
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            y = rand(1:teststate.bits)
            @fact swap(swap(teststate, x, y), y, x) --> roughly(teststate)
        end
    end
end

facts("CSWAP gate") do
    context("Three qubit mapping") do
        s000 = fromstates(QUBIT0, QUBIT0, QUBIT0)
        s001 = fromstates(QUBIT0, QUBIT0, QUBIT1)
        s010 = fromstates(QUBIT0, QUBIT1, QUBIT0)
        s011 = fromstates(QUBIT0, QUBIT1, QUBIT1)

        s100 = fromstates(QUBIT1, QUBIT0, QUBIT0)
        s101 = fromstates(QUBIT1, QUBIT0, QUBIT1)
        s110 = fromstates(QUBIT1, QUBIT1, QUBIT0)
        s111 = fromstates(QUBIT1, QUBIT1, QUBIT1)

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
        randomtest(3, 10) do teststate
            x = rand(1:teststate.bits - 2)
            y = rand(x + 1:teststate.bits - 1)
            z = rand(y + 1:teststate.bits)
            @fact cswap(cswap(teststate, x, y, z), x, z, y) --> roughly(teststate)
        end
    end
end

facts("Square root swap gate") do
    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
        @fact sqrtswap(s00, 1, 2) --> roughly(s00)
        @fact sqrtswap(s00, 2, 1) --> roughly(s00)

        @fact sqrtswap(s01, 1, 2) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0])
        @fact sqrtswap(s01, 2, 1) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0])

        @fact sqrtswap(s10, 1, 2) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0])
        @fact sqrtswap(s10, 2, 1) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0])

        @fact sqrtswap(s11, 1, 2) --> roughly(s11)
        @fact sqrtswap(s11, 2, 1) --> roughly(s11)
    end

    context("Three qubit mapping") do
        s000 = fromstates(QUBIT0, QUBIT0, QUBIT0)
        s001 = fromstates(QUBIT0, QUBIT0, QUBIT1)
        s010 = fromstates(QUBIT0, QUBIT1, QUBIT0)
        s011 = fromstates(QUBIT0, QUBIT1, QUBIT1)

        s100 = fromstates(QUBIT1, QUBIT0, QUBIT0)
        s101 = fromstates(QUBIT1, QUBIT0, QUBIT1)
        s110 = fromstates(QUBIT1, QUBIT1, QUBIT0)
        s111 = fromstates(QUBIT1, QUBIT1, QUBIT1)
        @fact sqrtswap(s000, 1, 2) --> roughly(s000)
        @fact sqrtswap(s000, 2, 1) --> roughly(s000)
        @fact sqrtswap(s000, 1, 3) --> roughly(s000)
        @fact sqrtswap(s000, 3, 1) --> roughly(s000)
        @fact sqrtswap(s000, 2, 3) --> roughly(s000)
        @fact sqrtswap(s000, 3, 2) --> roughly(s000)

        @fact sqrtswap(s001, 1, 2) --> roughly(s001)
        @fact sqrtswap(s001, 2, 1) --> roughly(s001)
        @fact sqrtswap(s001, 1, 3) --> roughly([0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrtswap(s001, 3, 1) --> roughly([0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrtswap(s001, 2, 3) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0, 0, 0, 0, 0])
        @fact sqrtswap(s001, 3, 2) --> roughly([0, 0.5 + 0.5im, 0.5 - 0.5im, 0, 0, 0, 0, 0])

        @fact sqrtswap(s010, 1, 2) --> roughly([0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrtswap(s010, 2, 1) --> roughly([0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0, 0])
        @fact sqrtswap(s010, 1, 3) --> roughly(s010)
        @fact sqrtswap(s010, 3, 1) --> roughly(s010)
        @fact sqrtswap(s010, 2, 3) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0, 0, 0, 0, 0])
        @fact sqrtswap(s010, 3, 2) --> roughly([0, 0.5 - 0.5im, 0.5 + 0.5im, 0, 0, 0, 0, 0])

        @fact sqrtswap(s011, 1, 2) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0])
        @fact sqrtswap(s011, 2, 1) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0.5 - 0.5im, 0, 0])
        @fact sqrtswap(s011, 1, 3) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0])
        @fact sqrtswap(s011, 3, 1) --> roughly([0, 0, 0, 0.5 + 0.5im, 0, 0, 0.5 - 0.5im, 0])
        @fact sqrtswap(s011, 2, 3) --> roughly(s011)
        @fact sqrtswap(s011, 3, 2) --> roughly(s011)

        @fact sqrtswap(s100, 1, 2) --> roughly([0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrtswap(s100, 2, 1) --> roughly([0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrtswap(s100, 1, 3) --> roughly([0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrtswap(s100, 3, 1) --> roughly([0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0, 0, 0])
        @fact sqrtswap(s100, 2, 3) --> roughly(s100)
        @fact sqrtswap(s100, 3, 2) --> roughly(s100)

        @fact sqrtswap(s101, 1, 2) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0])
        @fact sqrtswap(s101, 2, 1) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0.5 + 0.5im, 0, 0])
        @fact sqrtswap(s101, 1, 3) --> roughly(s101)
        @fact sqrtswap(s101, 3, 1) --> roughly(s101)
        @fact sqrtswap(s101, 2, 3) --> roughly([0, 0, 0, 0, 0, 0.5 + 0.5im, 0.5 - 0.5im, 0])
        @fact sqrtswap(s101, 3, 2) --> roughly([0, 0, 0, 0, 0, 0.5 + 0.5im, 0.5 - 0.5im, 0])

        @fact sqrtswap(s110, 1, 2) --> roughly(s110)
        @fact sqrtswap(s110, 2, 1) --> roughly(s110)
        @fact sqrtswap(s110, 1, 3) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0])
        @fact sqrtswap(s110, 3, 1) --> roughly([0, 0, 0, 0.5 - 0.5im, 0, 0, 0.5 + 0.5im, 0])
        @fact sqrtswap(s110, 2, 3) --> roughly([0, 0, 0, 0, 0, 0.5 - 0.5im, 0.5 + 0.5im, 0])
        @fact sqrtswap(s110, 3, 2) --> roughly([0, 0, 0, 0, 0, 0.5 - 0.5im, 0.5 + 0.5im, 0])

        @fact sqrtswap(s111, 1, 2) --> roughly(s111)
        @fact sqrtswap(s111, 2, 1) --> roughly(s111)
        @fact sqrtswap(s111, 1, 3) --> roughly(s111)
        @fact sqrtswap(s111, 3, 1) --> roughly(s111)
        @fact sqrtswap(s111, 2, 3) --> roughly(s111)
        @fact sqrtswap(s111, 3, 2) --> roughly(s111)
    end
end

facts("Phase shift gate") do
    context("Single qubit mapping") do
        @fact phaseshift(QUBIT0, 1, pi / 8) --> roughly(QUBIT0)
        @fact phaseshift(QUBIT1, 1, pi / 8) --> roughly([0.0 + 0im, exp(im * pi / 8)])
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)

        @fact phaseshift(s00, 1, pi / 8) --> roughly(s00)
        @fact phaseshift(s00, 2, pi / 8) --> roughly(s00)

        @fact phaseshift(s01, 1, pi / 8) --> roughly(s01)
        @fact phaseshift(s01, 2, pi / 8) --> roughly([0, exp(im * pi / 8), 0, 0])

        @fact phaseshift(s10, 1, pi / 8) --> roughly([0, 0, exp(im * pi / 8), 0])
        @fact phaseshift(s10, 2, pi / 8) --> roughly(s10)

        @fact phaseshift(s11, 1, pi / 8) --> roughly([0, 0, 0, exp(im * pi / 8)])
        @fact phaseshift(s11, 2, pi / 8) --> roughly([0, 0, 0, exp(im * pi / 8)])
    end

    context("Unitarity on N qubit states") do
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            theta = rand(1.0:pi)
            @fact phaseshift(phaseshift(teststate, x, theta), x, -theta) --> roughly(teststate)
        end
    end
end

facts("Pauli-X gate") do
    context("Single qubit mapping") do
        @fact paulix(QUBIT0, 1) --> roughly(QUBIT1)
        @fact paulix(QUBIT1, 1) --> roughly(QUBIT0)
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
        @fact paulix(s00, 1) --> roughly(s10)
        @fact paulix(s00, 2) --> roughly(s01)

        @fact paulix(s01, 1) --> roughly(s11)
        @fact paulix(s01, 2) --> roughly(s00)

        @fact paulix(s10, 1) --> roughly(s00)
        @fact paulix(s10, 2) --> roughly(s11)

        @fact paulix(s11, 1) --> roughly(s01)
        @fact paulix(s11, 2) --> roughly(s10)
    end

    context("Unitarity on N qubit states") do
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            @fact paulix(paulix(teststate, x), x) --> roughly(teststate)
        end
    end
end

facts("Pauli-Y gate") do
    context("Single qubit mapping") do
        @fact pauliy(QUBIT0, 1) --> roughly([0.0, -im])
        @fact pauliy(QUBIT1, 1) --> roughly([im, 0.0])
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
        @fact pauliy(s00, 1) --> roughly([0.0, 0, -im, 0])
        @fact pauliy(s00, 2) --> roughly([0.0, -im, 0, 0])

        @fact pauliy(s01, 1) --> roughly([0.0, 0, 0, -im])
        @fact pauliy(s01, 2) --> roughly([ im, 0, 0, 0.0])

        @fact pauliy(s10, 1) --> roughly([ im, 0, 0, 0.0])
        @fact pauliy(s10, 2) --> roughly([0.0, 0, 0, -im])

        @fact pauliy(s11, 1) --> roughly([0.0, im, 0, 0])
        @fact pauliy(s11, 2) --> roughly([0.0, 0, im, 0])
    end
end

facts("Pauli-Z gate") do
    context("Single qubit mapping") do
        @fact pauliz(QUBIT0, 1) --> roughly(QUBIT0)
        @fact pauliz(QUBIT1, 1) --> roughly([0.0 + 0im, -1])
    end

    context("Two qubit mapping") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)
        @fact pauliz(s00, 1) --> roughly(s00)
        @fact pauliz(s00, 2) --> roughly(s00)

        @fact pauliz(s01, 1) --> roughly(s01)
        @fact pauliz(s01, 2) --> roughly([0.0 + 0im, -1, 0, 0])

        @fact pauliz(s10, 1) --> roughly([0.0 + 0im, 0, -1, 0])
        @fact pauliz(s10, 2) --> roughly(s10)

        @fact pauliz(s11, 1) --> roughly([0.0 + 0im, 0, 0, -1])
        @fact pauliz(s11, 2) --> roughly([0.0 + 0im, 0, 0, -1])
    end

    context("Unitarity on N qubit states") do
        randomtest(2, 10) do teststate
            x = rand(1:teststate.bits)
            @fact pauliz(pauliz(teststate, x), x) --> roughly(teststate)
        end
    end
end

facts("Arbitrary unitary gates") do
    context("N-qubit Hadamard") do
        h = 1.0 / sqrt(2.0) * [1 1; 1 -1]
        randomtest(1, 5) do teststate
            x = rand(1:teststate.bits)
            @fact hadamard(teststate, x) --> roughly(unitary(teststate, h, x))
        end
    end

    context("N-qubit Not") do
        n = [0 1; 1 0]
        randomtest(1, 5) do teststate
            x = rand(1:teststate.bits)
            @fact Gates.not(teststate, x) --> roughly(unitary(teststate, n, x))
        end
    end

    context("N-qubit CNOT") do
        cn = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]
        randomtest(2, 5) do teststate
            x = rand(1:teststate.bits - 1)
            y = rand(x + 1:teststate.bits)
            @fact cnot(teststate, x, y) --> roughly(unitary(teststate, cn, x, y))
        end
    end

    context("N-qubit CCNOT") do
        ccn = [1 0 0 0 0 0 0 0;
               0 1 0 0 0 0 0 0;
               0 0 1 0 0 0 0 0;
               0 0 0 1 0 0 0 0;
               0 0 0 0 1 0 0 0;
               0 0 0 0 0 1 0 0;
               0 0 0 0 0 0 0 1;
               0 0 0 0 0 0 1 0]
        randomtest(3, 6) do teststate
            x = rand(1:teststate.bits - 2)
            y = rand(x + 1:teststate.bits - 1)
            z = rand(y + 1:teststate.bits)
            @fact ccnot(teststate, z, x, y) --> roughly(unitary(teststate, ccn, z, x, y))
        end
    end

    context("N-qubit SWAP") do
        s = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
        randomtest(2, 5) do teststate
            x = rand(1:teststate.bits - 1)
            y = rand(x + 1:teststate.bits)
            @fact swap(teststate, x, y) --> roughly(unitary(teststate, s, x, y))
        end
    end
    context("N-qubit CSWAP") do
        cs = [1 0 0 0 0 0 0 0;
              0 1 0 0 0 0 0 0;
              0 0 1 0 0 0 0 0;
              0 0 0 1 0 0 0 0;
              0 0 0 0 1 0 0 0;
              0 0 0 0 0 0 1 0;
              0 0 0 0 0 1 0 0;
              0 0 0 0 0 0 0 1]
        randomtest(3, 6) do teststate
            x = rand(1:teststate.bits - 2)
            y = rand(x + 1:teststate.bits - 1)
            z = rand(y + 1:teststate.bits)
            @fact cswap(teststate, x, y, z) --> roughly(unitary(teststate, cs, x, y, z))
        end
    end
end

function todecimal(darr::Vector{Int})
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
        @fact measure(QUBIT0) --> (EMPTY_STATE, [0])
        @fact measure(QUBIT1) --> (EMPTY_STATE, [1])
    end

    context("Equal probability statistical test") do
        s0 = fromvector([0.5, -0.5, 0 - 0.5im, 0 + 0.5im])

        nruns = 10000
        stats = [0, 0, 0, 0]
        for i = 1:nruns
            posterior, measurement = measure(s0)
            n = todecimal(measurement)
            @fact posterior[n + 1] --> 1.0 + 0.0im
            @fact countnz(posterior.vector) --> 1
            stats[n + 1] += 1
        end
        @fact stats[1] / nruns --> roughly(0.25, atol=0.1)
        @fact stats[2] / nruns --> roughly(0.25, atol=0.1)
        @fact stats[3] / nruns --> roughly(0.25, atol=0.1)
        @fact stats[4] / nruns --> roughly(0.25, atol=0.1)
    end

    context("Null probability statistical test") do
        rsqrt8 = 1.0 / sqrt(8)
        s0 = fromvector([rsqrt8 + rsqrt8 * im, rsqrt8 - rsqrt8 * im, 0.0,  1.0 / sqrt(2)])

        nruns = 10000
        stats = [0, 0, 0, 0]
        for i=1:nruns
            posterior, measurement = measure(s0)
            n = todecimal(measurement)
            @fact posterior[n + 1] --> 1.0 + 0.0im
            @fact countnz(posterior.vector) --> 1
            stats[n + 1] += 1
        end
        @fact stats[1] / nruns --> roughly(0.25, atol=0.1)
        @fact stats[2] / nruns --> roughly(0.25, atol=0.1)
        @fact stats[3] / nruns --> 0
        @fact stats[4] / nruns --> roughly(0.5, atol=0.1)
    end
end

# These tests are not that comprehensive, but they are hard to automate
facts("Partial Measurement") do
    context("Single qubit measurement") do
        @fact partialmeasure(QUBIT0, 1) --> (EMPTY_STATE, [0])
        @fact partialmeasure(QUBIT1, 1) --> (EMPTY_STATE, [1])
    end

    context("Two qubit measurement") do
        s00 = fromstates(QUBIT0, QUBIT0)
        s01 = fromstates(QUBIT0, QUBIT1)
        s10 = fromstates(QUBIT1, QUBIT0)
        s11 = fromstates(QUBIT1, QUBIT1)

        # State |00>
        posterior, measurement = partialmeasure(s00, 1)
        @fact measurement --> [0]
        @fact posterior --> roughly(QUBIT0)

        posterior, measurement = partialmeasure(s00, 2)
        @fact measurement --> [0]
        @fact posterior --> roughly(QUBIT0)

        posterior, measurement = partialmeasure(s00, 1, 2)
        @fact measurement --> [0, 0]
        @fact posterior --> roughly(EMPTY_STATE)

        # State |01>
        posterior, measurement = partialmeasure(s01, 1)
        @fact measurement --> [0]
        @fact posterior --> roughly(QUBIT1)

        posterior, measurement = partialmeasure(s01, 2)
        @fact measurement --> [1]
        @fact posterior --> roughly(QUBIT0)

        posterior, measurement = partialmeasure(s01, 1, 2)
        @fact measurement --> [0, 1]
        @fact posterior --> roughly(EMPTY_STATE)

        # State |10>
        posterior, measurement = partialmeasure(s10, 1)
        @fact measurement --> [1]
        @fact posterior --> roughly(QUBIT0)

        posterior, measurement = partialmeasure(s10, 2)
        @fact measurement --> [0]
        @fact posterior --> roughly(QUBIT1)

        posterior, measurement = partialmeasure(s10, 1, 2)
        @fact measurement --> [1, 0]
        @fact posterior --> roughly(EMPTY_STATE)

        # State |11>
        posterior, measurement = partialmeasure(s11, 1)
        @fact measurement --> [1]
        @fact posterior --> roughly(QUBIT1)

        posterior, measurement = partialmeasure(s11, 2)
        @fact measurement --> [1]
        @fact posterior --> roughly(QUBIT1)

        posterior, measurement = partialmeasure(s11, 1, 2)
        @fact measurement --> [1, 1]
        @fact posterior --> roughly(EMPTY_STATE)

        # Try states with different probabilities
        h = hadamard(s00, 1)
        for i=1:100
            posterior, measurement = partialmeasure(h, 1)
            @fact posterior --> roughly(QUBIT0)
            posterior, measurement = partialmeasure(h, 2)
            @fact measurement --> [0]
            @fact posterior --> roughly(1.0 / sqrt(2) .* [1.0 + 0im, 1])
        end

        # Partially measuring the same bit multiple time should have the same effect
        # as measuring it only once
        p1, m1 = partialmeasure(s11, 1, 1, 1)
        p2, m2 = partialmeasure(s11, 1)
        @fact p1 --> roughly(p2)
        @fact m2 --> m2

        p1, m1 = partialmeasure(s11, 1, 1, 2, 2, 1, 2)
        p2, m2 = partialmeasure(s11, 1, 2)
        @fact p1 --> roughly(p2)
        @fact m2 --> m2
    end
end
