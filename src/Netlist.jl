module Netlist

using ParserCombinator
using ..State, ..Gates

export Gate, InitialState, parse_netlist

const GATE_MAP = Dict{ASCIIString, Function}(
    "hadamard"    => hadamard,
    "not"         => not,
    "cnot"        => cnot,
    "ccnot"       => ccnot,
    "phase_shift" => phase_shift,
    "pauli_x"     => pauli_x,
    "pauli_y"     => pauli_y,
    "pauli_z"     => pauli_z,
    "swap"        => swap,
    "cswap"       => cswap,
    "sqrt_swap"   => sqrt_swap,
    "probe"       => probe,
    "measure"     => measure,
    "partial_measure" => partial_measure
)

type Gate
    index::Int
    fn::Function
    inputs::Vector{Int}
    arguments::Vector{Any}
end

type InitialState
    index::Int
    state::QuantumState
end

function parse_netlist{S<:AbstractString}(filename::S)
    file = open(filename)
    netlist_string = readall(file)
    close(file)

    spc = Drop(Star(Space()))

    @with_pre spc begin
        num = PInt64() | PFloat64()
        str = p"\"[\w ]*\""
        con = e"pi" | e"i" > parse_constant
        op  = e"+" | e"-" | e"*" | e"/"

        multiplied_constant = num + con > multiply_constant
        numeric = multiplied_constant | con | num
        simple_expr = (numeric + op + numeric) | numeric > perform_op

        index = PInt64() + spc

        # First the grammar describing the gates
        gate_name  = e"hadamard"  | e"not"       | e"cnot"            | e"ccnot"   | e"swap"    |
                     e"cswap"     | e"sqrt_swap" | e"phase_shift"     | e"pauli_x" | e"pauli_y" |
                     e"pauli_z"   | e"measure"   | e"partial_measure" | e"probe" > parse_gate_name

        conn_list = E"[]" | (E"[" + index + Star(E"," + index) + E"]") > connections

        args      = (str | simple_expr | numeric) + spc
        args_list = E"()" | (E"(" + args + Star(E"," + args) + E")") > arguments

        gate = index + E":" + gate_name + ((conn_list + args_list) | (conn_list | args_list)) > create_gate

        # Then the grammar for the initial quantum state
        state_name = E"qstate"
        qubit = (e"bell" | PInt64()) + spc > parse_qubit
        qubit_list = E"(" + qubit + Star(E"," + qubit) + E")" > from_states
        initial_state = index + E":" + state_name + qubit_list > InitialState

        all_gates = Plus(gate | initial_state) + spc + Eos()
    end

    netlist = parse_one(netlist_string, all_gates)

    gates::Vector{Gate}            = filter(x -> isa(x, Gate), netlist)
    states::Vector{InitialState}   = filter(x -> isa(x, InitialState), netlist)
    outputs::Vector{Nullable{Any}} = fill(Nullable{Any}(), length(netlist))
    return states, gates, outputs
end

function parse_constant{S<:AbstractString}(constant::S)
    if constant == "pi"
        return pi
    elseif constant == "i"
        return im
    end
    error("Unsupported constant in expression")
end

function multiply_constant(multiplier, constant)
    return multiplier * constant
end

function perform_op{S<:AbstractString}(num1, op::S, num2)
    if op == "+"
        return num1 + num2
    elseif op == "-"
        return num1 - num2
    elseif op == "*"
        return num1 * num2
    elseif op == "/"
        return num1 / num2
    end
    error("Unsupported operator in expresseion")
end

perform_op(num) = num

function parse_gate_name{S<:AbstractString}(gate::S)
    if haskey(GATE_MAP, gate)
       return GATE_MAP[gate]
    end
    error("Unsupported gate ($gate) in netlist")
end

function parse_qubit{S<:AbstractString}(qubit::S)
    if qubit == "bell"
        return BELL_STATE
    end
    error("Unsupported qubit state in arguments")
end

function parse_qubit(qubit::Int)
    if qubit == 0
        return QUBIT_0
    elseif qubit == 1
        return QUBIT_1
    end
    error("Unsupported qubit state in arguments")
end

function connections(cs::Int...)
    return Vector{Int}([cs...])
end

function arguments(args::Any...)
    return Vector{Any}([args...])
end

function arguments()
    return Vector{Any}([])
end

function create_gate(index::Int, fn::Function, connections::Vector{Int})
    return Gate(index, fn, connections, Vector{Any}([]))
end

create_gate(index::Int, fn::Function, connections::Vector{Int}, args::Vector{Any}) = Gate(index, fn, connections, args)


end
