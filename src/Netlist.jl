module Netlist

using Iterators

using QSpice.State, QSpice.Gates

export Gate, parse_netlist

const NAME_FN_MAP = Dict(
    "hadamard"    => hadamard,
    "cnot"        => cnot,
    "not"         => not,
    "phase_shift" => phase_shift,
    "pauli_x"     => pauli_x,
    "pauli_y"     => pauli_y,
    "pauli_z"     => pauli_z,
    "swap"        => swap,
    "probe"       => probe,
    "measure"     => measure,
    "partial_measure" => partial_measure
)

type Gate
    index::Int
    fn::Function
    inputs::Vector{Int}
    arguments::Vector{Int}
end

type InitialState
    index::Int
    state::QuantumState
end

function read_gates(filename::AbstractString)
    file = open(filename);
    gate_strings = imap(l -> split(l), eachline(file)) |> collect
    close(file)
    return vcat(gate_strings...)
end

function parse_one(gate_string::AbstractString)
    index, rest = split(gate_string, ":")
    index = parse(Int, index)

    name, rest = split(rest, "[")

    if name == "qstate"
        args = split(split(rest, "(")[2], ")")[1]
        args = split(args, ",")
        return InitialState(index, from_states(map(x -> x == "0" ? QUBIT_0 : QUBIT_1, args)...))
        return gate_string
    end

    fn = NAME_FN_MAP[name]

    inputs, rest = split(rest, "]")
    inputs = map(x -> parse(Int, x), split(inputs, ","))

    args = split(split(rest, "(")[2], ")")[1]
    if length(args) > 0
        args = map(x -> parse(Int, x), split(args, ","))
    else
        args = []
    end
    return Gate(index, fn, inputs, args)
end

function parse_netlist(filename::AbstractString)
    gate_strings = read_gates(filename)
    outputs = fill(Nullable{QuantumState}(), length(gate_strings))
    elements = map(parse_one, gate_strings)

    initial_states = filter(x -> isa(x, InitialState), elements)
    filter!(x -> isa(x, Gate), elements)

    gates::Vector{Gate} = elements
    return (initial_states, gates, outputs)
end

end
