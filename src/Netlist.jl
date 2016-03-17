module Netlist

using QSpice.State, QSpice.Gates

export Gate, parse_netlist

const NAME_FN_MAP = Dict{ASCIIString, Function}(
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

function ignore_whitespace{S<:AbstractString}(stream::S)
    if stream == ""
        return stream
    end

    pos = 1
    while isspace(stream[pos]) pos += 1 end
    return stream[pos:end]
end

function parse_index{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream =="" || !isdigit(stream[1])
        return (Nullable{Int}(), stream)
    end

    index_end = findfirst(stream, ':')
    return (tryparse(Int, stream[1:index_end - 1]), stream[index_end + 1:end])
end

function parse_name{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream == "" || !isalpha(stream[1])
        return (Nullable{AbstractString}(), stream)
    end

    bracket = findfirst(stream, '[')
    paren = findfirst(stream, '(')

    if bracket == 0 && paren == 0
        return (Nullable{AbstractString}(), stream)
    end

    index_end = bracket == 0 ? paren : (paren == 0 ? bracket : min(paren, bracket))
    return (Nullable(stream[1:index_end - 1]), stream[index_end:end])
end

function parse_connections{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream == "" || stream[1] != '['
        return (Nullable{Array{Int}}(), stream)
    end

    index_end = findfirst(stream, ']')
    split_connections = split(stream[2:index_end - 1], ',', keep = false)

    if isempty(split_connections)
        return (Nullable{Array{Int}}(), stream[index_end + 1:end])
    end

    connections = map(x -> parse(Int, strip(x)), split_connections)
    return (Nullable(connections), stream[index_end + 1:end])
end

function tryparse_floatint{S<:AbstractString}(arg::S)
    # First try a simple integer parser
    maybe_int = tryparse(Int, arg)
    if !isnull(maybe_int)
        return Nullable(get(maybe_int))
    end

    # If that didn't work, try a simple float parser
    return tryparse(Float64, arg)
end

function tryparse_special{S<:AbstractString, T<:AbstractString, U<:Number}(arg::S, special::T, value::U)
    sarg = split(arg, special, keep = false)
    if isempty(sarg)
        return Nullable(value)
    end

    if length(sarg) != 1
        return Nullable{Number}()
    end

    maybe_numeric = tryparse_floatint(sarg[1])
    return isnull(maybe_numeric) ? Nullable{Number}() : Nullable(get(maybe_numeric) * value)
end

# A very primitive expression parser. It supports *one of* +-*/,
# the special constants i and pi, floating point and integer numbers
function tryparse_expression{S<:AbstractString}(arg::S)
    op = max(findfirst(arg, '+'),
             findfirst(arg, '-'),
             findfirst(arg, '*'),
             findfirst(arg, '/'))

    if op == 0
        return Nullable{Number}()
    end

    # Now we only support a single operator
    sargs = split(arg, arg[op])
    if length(sargs) != 2
        return Nullable{Number}()
    end

    numargs = []
    for a in sargs
        a = strip(a)
        result = tryparse_floatint(a)
        if contains(a, "pi")
            result = tryparse_special(a, "pi", pi)
        elseif contains(a, "i")
            result = tryparse_special(a, "i", im)
        end

        if isnull(result)
            return Nullable{Number}()
        end

        push!(numargs, get(result))
    end

    if arg[op] == '+'
        return Nullable(numargs[1] + numargs[2])
    elseif arg[op] == '-'
        return Nullable(numargs[1] - numargs[2])
    elseif arg[op] == '*'
        return Nullable(numargs[1] * numargs[2])
    elseif arg[op] == '/'
        return Nullable(numargs[1] / numargs[2])
    end

    return Nullable{Number}()
end

function parse_argument{S<:AbstractString}(arg::S)
    arg = strip(arg)

    maybe_numeric = tryparse_floatint(arg)
    if !isnull(maybe_numeric)
        return get(maybe_numeric)
    end

    # Now try and parse it as a supported expression
    maybe_expr = tryparse_expression(arg)

    # If nothing worked, we'll treat it as a string
    return isnull(maybe_expr) ? arg : get(maybe_expr)
end

function parse_arguments{S<:AbstractString}(stream::S)
    stream = ignore_whitespace(stream)

    if stream == "" || stream[1] != '('
        return (Nullable{Array{Any}}(), stream)
    end

    index_end = findfirst(stream, ')')
    split_args = split(stream[2:index_end - 1], ',', keep = false)

    if isempty(split_args)
        return (Nullable{Array{Any}}(), stream[index_end + 1: end])
    end

    arguments::Array{Any} = map(parse_argument, split_args)
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

            actual_args = isnull(arguments) ? [] : get(arguments)
            push!(gates, Gate(get(index), fn, get(connections), actual_args))
        end
    end
    outputs = fill(Nullable{Any}(), length(initial_state) + length(gates))
    return (initial_state, gates, outputs)
end

end
