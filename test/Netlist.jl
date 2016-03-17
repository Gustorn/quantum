include("../src/QSpice.jl")
include("../src/Netlist.jl")
include("Common.jl")

using FactCheck

using Netlist, QSpice.State, QSpice.Gates

function result_equals{T, U}(t1::Tuple{Nullable{T}, U})
    function detail(t2::Tuple{Nullable{T}, U})
        if isnull(t1[1]) && isnull(t2[1])
            return t1[2] == t2[2]
        elseif isnull(t1[1]) != isnull(t2[1])
            return false
        else
            return get(t1[1]) == get(t2[1]) && t1[2] == t2[2]
        end
    end
    return detail
end

facts("Test individual parsing functions") do
    context("Ignoring whitespace") do
        @fact Netlist.ignore_whitespace("") --> ""
        @fact Netlist.ignore_whitespace("test") --> "test"
        @fact Netlist.ignore_whitespace("   test") --> "test"
        @fact Netlist.ignore_whitespace("test   ") --> "test   "
        @fact Netlist.ignore_whitespace("   test   ") --> "test   "
        @fact Netlist.ignore_whitespace("   test test   test   ") --> "test test   test   "
    end

    context("Parsing index") do
        @fact Netlist.parse_index("") --> result_equals((Nullable{Int}(), ""))
        @fact Netlist.parse_index("5:") --> result_equals((Nullable(5), ""))
        @fact Netlist.parse_index("   10:") --> result_equals((Nullable(10), ""))
        @fact Netlist.parse_index("100:hadamard") --> result_equals((Nullable(100), "hadamard"))
        @fact Netlist.parse_index("h10:hadamard") --> result_equals((Nullable{Int}(), "h10:hadamard"))
    end

    context("Parsing gate name") do
        @fact Netlist.parse_name("") --> result_equals((Nullable{AbstractString}(), ""))
        @fact Netlist.parse_name("hadamard") --> result_equals((Nullable{AbstractString}(), "hadamard"))
        @fact Netlist.parse_name("hadamard[10") --> result_equals((Nullable("hadamard"), "[10"))
        @fact Netlist.parse_name("cnot(10") --> result_equals((Nullable("cnot"), "(10"))
        @fact Netlist.parse_name("swap[10](20)") --> result_equals((Nullable("swap"), "[10](20)"))
        @fact Netlist.parse_name("  hadamard[10") --> result_equals((Nullable("hadamard"), "[10"))
    end

    context("Parsing connections") do
        @fact Netlist.parse_connections("") --> result_equals((Nullable{Array{Int}}(), ""))
        @fact Netlist.parse_connections("[]") --> result_equals((Nullable{Array{Int}}(), ""))
        @fact Netlist.parse_connections("[1]") --> result_equals((Nullable([1]), ""))
        @fact Netlist.parse_connections("[1](2,3)") --> result_equals((Nullable([1]), "(2,3)"))
        @fact Netlist.parse_connections("[1,2](3,4)") --> result_equals((Nullable([1,2]), "(3,4)"))
        @fact Netlist.parse_connections("[1,2](Test [])") --> result_equals((Nullable([1,2]), "(Test [])"))
        @fact Netlist.parse_connections("[   1,  2   ](3,4)") --> result_equals((Nullable([1, 2]), "(3,4)"))
        @fact Netlist.parse_connections("  [1,2](3,4) 10:hadamard") --> result_equals((Nullable([1, 2]), "(3,4) 10:hadamard"))
    end

    context("Parsing arguments") do
        @fact Netlist.parse_arguments("") --> result_equals((Nullable{Array{Any}}(), ""))
        @fact Netlist.parse_arguments("()") --> result_equals((Nullable{Array{Any}}(), ""))

        # Testing simple arguments
        expected::Vector{Any} = [10]
        @fact Netlist.parse_arguments("(10)") --> result_equals((Nullable(expected), ""))

        expected = ["Test: String[]?!"]
        @fact Netlist.parse_arguments("(Test: String[]?!)") --> result_equals((Nullable(expected), ""))
        @fact Netlist.parse_arguments("(10, Test String, 20) 10:swap") --> result_equals((Nullable([10, "Test String", 20]), " 10:swap"))
        @fact Netlist.parse_arguments("   (  10, 20, Foo  )  a") --> result_equals((Nullable([10, 20, "Foo"]), "  a"))
        @fact Netlist.parse_arguments("(11.46, 10, Foo)") --> result_equals((Nullable([11.46, 10, "Foo"]), ""))


        # Testing expressions
        @fact Netlist.parse_arguments("(10 + 20, T) a") --> result_equals((Nullable([30, "T"]), " a"))
        @fact Netlist.parse_arguments("(10 + 20 + 30, 10)") --> result_equals((Nullable(["10 + 20 + 30", 10]), ""))
        @fact Netlist.parse_arguments("(pi/8, T)") --> result_equals((Nullable([pi/8, "T"]), ""))
        @fact Netlist.parse_arguments("(0.5pi/8, T)") --> result_equals((Nullable([0.5 * pi / 8, "T"]), ""))
        @fact Netlist.parse_arguments("(5pi - 10i, T)") --> result_equals((Nullable([5 * pi - 10im, "T"]), ""))
        @fact Netlist.parse_arguments("(5i * 5.0i, T)") --> result_equals((Nullable([5im * 5.0im, "T"]), ""))
    end

    context("Qubit parsing") do
        @fact Netlist.parse_qubit(0) --> QUBIT_0
        @fact Netlist.parse_qubit(1) --> QUBIT_1
        @fact Netlist.parse_qubit("bell") --> BELL_STATE
        @fact Netlist.parse_qubit("BELL") --> BELL_STATE
    end
end

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
        check_gate(13, Gates.probe,           [12], Vector{Any}(["Foo"]))
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
