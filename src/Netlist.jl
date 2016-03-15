module Netlist

using Iterators

using QSpice.State, QSpice.Gates

export Gate, parse_netlist

const NAME_FN_MAP = Dict{ASCIIString, Function}(
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
    arguments::Vector{Any}
end

type InitialState
    index::Int
    state::QuantumState
end

function ignore_whitespace{S<:AbstractString}(stream::S)
    pos = 1
    while isspace(stream[pos]) pos += 1 end
    return stream[pos:end]
end

function parse_index{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if !isdigit(stream[1])
        return (Nullable{Int}(), stream)
    end

    index_end = findfirst(stream, ':')
    return (tryparse(Int, stream[1:index_end - 1]), stream[index_end + 1:end])
end

function parse_name{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if !isalpha(stream[1])
        return (Nullable{AbstractString}(), stream)
    end

    index_end = min(findfirst(stream, '['), findfirst(stream, '('))
    return (Nullable(stream[1:index_end - 1]), stream[index_end:end])
end

function parse_connections{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream[1] != '['
        return (Nullable{Array{Int}}(), stream)
    end

    index_end = findfirst(stream, ']')
    split_connections = split(stream[2:index_end - 1], ',', keep = false)
    connections = imap(x -> parse(Int, strip(x)), split_connections) |> collect
    return (Nullable(connections), stream[index_end + 1:end])
end

function parse_argument{S<:AbstractString}(arg::S)
    arg = strip(arg)
    maybe_int = tryparse(Int, arg)
    return isnull(maybe_int) ? arg : get(maybe_int)
end

function parse_arguments{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream[1] != '('
        return (Nullable{Array{Any}}(), stream)
    end

    index_end = findfirst(stream, ')')
    split_args = split(stream[2:index_end - 1], ',', keep = false)
    arguments = imap(parse_argument, split_args) |> collect
    return (Nullable(arguments), stream[index_end + 1:end])
end

function parse_qubit{S<:AbstractString}(qubit::S)
    if lowercase(qubit) == "bell"
        return BELL_STATE
    end
    error("Unsupported qubit state")
end

function parse_qubit(qubit::Int)
    if qubit == 0
        return QUBIT_0
    elseif qubit == 1
        return QUBIT_1
    end
    error("Unsupported qubit state")
end

function parse_netlist{S<:AbstractString}(filename::S)
    file = open(filename)
    netlist = strip(readall(file))
    close(file)

    initial_state::Vector{InitialState} = []
    gates::Vector{Gate} = []

    while netlist != ""
        index, netlist = parse_index(netlist)
        if isnull(index) error("No index found for gate in the netlist") end

        name, netlist = parse_name(netlist)
        if isnull(name) error("No name found for gate in the netlist") end

        connections, netlist = parse_connections(netlist)
        arguments, netlist = parse_arguments(netlist)

        if get(name) == "qstate"
            if isnull(arguments) error("No qubits given for an initial quantum state") end
            push!(initial_state, InitialState(get(index), from_states(map(parse_qubit, get(arguments))...)))
        else
            if !haskey(NAME_FN_MAP, get(name))
                error("The netlist contains a quantum gate that's not implemented yet")
            end
            fn = NAME_FN_MAP[get(name)]
            push!(gates, Gate(get(index), fn, get(connections), get(arguments)))
        end
    end
    outputs = fill(Nullable{QuantumState}(), length(initial_state) + length(gates))
    return (initial_state, gates, outputs)
end

end
