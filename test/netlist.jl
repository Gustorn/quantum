include("common.jl")

using FactCheck

using Netlist, QSpice.State, QSpice.Gates

function resultequals{T, U}(t1::Tuple{T, Nullable{U}})
    function detail(t2::Tuple{T, Nullable{U}})
        if isnull(t1[2]) && isnull(t2[2])
            return t1[1] == t2[1]
        elseif isnull(t1[2]) != isnull(t2[2])
            return false
        else
            return get(t1[2]) == get(t2[2]) && t1[1] == t2[1]
        end
    end
    return detail
end

facts("Test individual parsing functions") do
    context("Ignoring whitespace") do
        @fact Netlist.skipwhitespace("") --> ""
        @fact Netlist.skipwhitespace("test") --> "test"
        @fact Netlist.skipwhitespace("   test") --> "test"
        @fact Netlist.skipwhitespace("test   ") --> "test   "
        @fact Netlist.skipwhitespace("   test   ") --> "test   "
        @fact Netlist.skipwhitespace("   test test   test   ") --> "test test   test   "
    end

    context("Parsing index") do
        @fact Netlist.gateoutput("") --> resultequals(("", Nullable{Int}()))
        @fact Netlist.gateoutput("5:") --> resultequals(("", Nullable(5)))
        @fact Netlist.gateoutput("   10:") --> resultequals(("", Nullable(10)))
        @fact Netlist.gateoutput("100:hadamard") --> resultequals(("hadamard", Nullable(100)))
        @fact Netlist.gateoutput("h10:hadamard") --> resultequals(("h10:hadamard", Nullable{Int}()))
    end

    context("Parsing gate name") do
        @fact Netlist.gatename("") --> resultequals(("", Nullable{ASCIIString}()))
        @fact Netlist.gatename("hadamard") --> resultequals(("hadamard", Nullable{ASCIIString}()))
        @fact Netlist.gatename("hadamard[10") --> resultequals(("[10", Nullable("hadamard")))
        @fact Netlist.gatename("cnot(10") --> resultequals(("(10", Nullable("cnot")))
        @fact Netlist.gatename("swap[10](20)") --> resultequals(("[10](20)", Nullable("swap")))
        @fact Netlist.gatename("  hadamard[10") --> resultequals(("[10", Nullable("hadamard")))
    end

    context("Parsing connections") do
        @fact Netlist.inputs("") --> resultequals(("", Nullable{Vector{Int}}()))
        @fact Netlist.inputs("[]") --> resultequals(("", Nullable{Vector{Int}}()))
        @fact Netlist.inputs("[1]") --> resultequals(("", Nullable([1])))
        @fact Netlist.inputs("[1](2,3)") --> resultequals(("(2,3)", Nullable([1])))
        @fact Netlist.inputs("[1,2](3,4)") --> resultequals(("(3,4)", Nullable([1,2])))
        @fact Netlist.inputs("[1,2](Test [])") --> resultequals(("(Test [])", Nullable([1,2])))
        @fact Netlist.inputs("[   1,  2   ](3,4)") --> resultequals(("(3,4)", Nullable([1, 2])))
        @fact Netlist.inputs("  [1,2](3,4) 10:hadamard") --> resultequals(("(3,4) 10:hadamard", Nullable([1, 2])))
    end

    context("Parsing arguments") do
        @fact Netlist.arguments("") --> resultequals(("", Nullable{Vector{Any}}()))
        @fact Netlist.arguments("()") --> resultequals(("", Nullable{Vector{Any}}()))

        # Testing simple arguments
        expected::Vector{Any} = [10]
        @fact Netlist.arguments("(10)") --> resultequals(("", Nullable(expected)))

        expected = ["Test: String[]?!"]
        @fact Netlist.arguments("(Test: String[]?!)") --> resultequals(("", Nullable(expected)))
        @fact Netlist.arguments("(10, Test String, 20) 10:swap") --> resultequals((" 10:swap", Nullable([10, "Test String", 20])))
        @fact Netlist.arguments("   (  10, 20, Foo  )  a") --> resultequals(("  a", Nullable([10, 20, "Foo"])))
        @fact Netlist.arguments("(11.46, 10, Foo)") --> resultequals(("", Nullable([11.46, 10, "Foo"])))


        # Testing expressions
        @fact Netlist.arguments("(10 + 20, T) a") --> resultequals((" a", Nullable([30, "T"])))
        @fact Netlist.arguments("(10 + 20 + 30, 10)") --> resultequals(("", Nullable(["10 + 20 + 30", 10])))
        @fact Netlist.arguments("(pi/8, T)") --> resultequals(("", Nullable([pi/8, "T"])))
        @fact Netlist.arguments("(0.5pi/8, T)") --> resultequals(("", Nullable([0.5 * pi / 8, "T"])))
        @fact Netlist.arguments("(5pi - 10i, T)") --> resultequals(("", Nullable([5 * pi - 10im, "T"])))
        @fact Netlist.arguments("(5i * 5.0i, T)") --> resultequals(("", Nullable([5im * 5.0im, "T"])))
    end

    context("Qubit parsing") do
        @fact Netlist.qubit(0) --> QUBIT0
        @fact Netlist.qubit(1) --> QUBIT1
        @fact Netlist.qubit("bell") --> BELL_STATE
        @fact Netlist.qubit("BELL") --> BELL_STATE
    end
end

facts("Testing a full example") do
    initial, gates, outputs = Netlist.parsenetlist("test/test.qnl")

    context("Checking if the initial state is properly initialized") do
        @fact length(initial) --> 1
        @fact initial[1].output --> 1
        @fact initial[1].state --> roughly(fromstates(BELL_STATE, QUBIT1, QUBIT1, QUBIT0))
    end

    context("Checking if all gates are properly read") do
        @fact length(gates) --> 16
        @fact length(unique(map(g -> g.output, gates))) --> length(gates)

        function checkgate(index::Int, fn::Function, inputs::Vector{Int}, arguments::Vector{Any})
            i = findfirst(x -> x.output == index, gates)
            @fact i --> FactCheck.not(0)
            @fact gates[i].fn --> exactly(fn)
            @fact gates[i].inputs --> inputs
            @fact gates[i].arguments --> arguments
        end

        checkgate(2,  Gates.hadamard,       [1],  Vector{Any}([1]))
        checkgate(3,  Gates.not,            [2],  Vector{Any}([1]))
        checkgate(4,  Gates.cnot,           [3],  Vector{Any}([1, 2]))
        checkgate(5,  Gates.ccnot,          [4],  Vector{Any}([1, 2, 3]))
        checkgate(6,  Gates.phaseshift,     [5],  Vector{Any}([1, pi/4]))
        checkgate(7,  Gates.paulix,         [6],  Vector{Any}([1]))
        checkgate(8,  Gates.pauliy,         [7],  Vector{Any}([1]))
        checkgate(9,  Gates.pauliz,         [8],  Vector{Any}([1]))
        checkgate(10, Gates.swap,           [9],  Vector{Any}([1,2]))
        checkgate(11, Gates.cswap,          [10], Vector{Any}([1,2,3]))
        checkgate(12, Gates.sqrtswap,       [11], Vector{Any}([2,3]))
        checkgate(13, Gates.probe,          [12], Vector{Any}(["Foo"]))
        checkgate(14, Gates.probe,          [12], Vector{Any}([]))
        checkgate(15, Gates.measure,        [12], Vector{Any}([]))
        checkgate(16, Gates.partialmeasure, [12], Vector{Any}([1]))
        checkgate(17, Gates.partialmeasure, [12], Vector{Any}([1,2,3]))
    end

    context("Checking if we have the correct number of outputs") do
        @fact length(outputs) --> length(gates) + length(initial)
        @fact all(isnull, outputs) --> true
    end
end
