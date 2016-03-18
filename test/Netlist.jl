include("Common.jl")

using FactCheck

using Netlist, QSpice.State, QSpice.Gates

facts("Testing a full example") do
    initial, gates, outputs = Netlist.parse_netlist("test/test.qnl")

    context("Checking if the initial state is properly initialized") do
        @fact length(initial) --> 1
        @fact initial[1].index --> 1
        @fact initial[1].state --> roughly(from_states(BELL_STATE, QUBIT_1, QUBIT_1, QUBIT_0))
    end

    context("Checking if all gates are properly read") do
        @fact length(gates) --> 16
        @fact length(unique(map(g -> g.index, gates))) --> length(gates)

        function check_gate(index::Int, fn::Function, inputs::Vector{Int}, arguments::Vector{Any})
            i = findfirst(x -> x.index == index, gates)
            @fact i --> FactCheck.not(0)
            @fact gates[i].fn --> exactly(fn)
            @fact gates[i].inputs --> inputs
            @fact gates[i].arguments --> arguments
        end

        check_gate(2,  Gates.hadamard,        [1],  Vector{Any}([1]))
        check_gate(3,  Gates.not,             [2],  Vector{Any}([1]))
        check_gate(4,  Gates.cnot,            [3],  Vector{Any}([1, 2]))
        check_gate(5,  Gates.ccnot,           [4],  Vector{Any}([1, 2, 3]))
        check_gate(6,  Gates.phase_shift,     [5],  Vector{Any}([1, pi/4]))
        check_gate(7,  Gates.pauli_x,         [6],  Vector{Any}([1]))
        check_gate(8,  Gates.pauli_y,         [7],  Vector{Any}([1]))
        check_gate(9,  Gates.pauli_z,         [8],  Vector{Any}([1]))
        check_gate(10, Gates.swap,            [9],  Vector{Any}([1,2]))
        check_gate(11, Gates.cswap,           [10], Vector{Any}([1,2,3]))
        check_gate(12, Gates.sqrt_swap,       [11], Vector{Any}([2,3]))
        check_gate(13, Gates.probe,           [12], Vector{Any}(["\"Foo\""]))
        check_gate(14, Gates.probe,           [12], Vector{Any}([]))
        check_gate(15, Gates.measure,         [12], Vector{Any}([]))
        check_gate(16, Gates.partial_measure, [12], Vector{Any}([1]))
        check_gate(17, Gates.partial_measure, [12], Vector{Any}([1,2,3]))
    end

    context("Checking if we have the correct number of outputs") do
        @fact length(outputs) --> length(gates) + length(initial)
        @fact all(isnull, outputs) --> true
    end
end
