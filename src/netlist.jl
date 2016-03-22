module Netlist

using QSpice.State, QSpice.Gates

export Gate, parsenetlist

const FUNCTION_MAP = Dict{ASCIIString, Function}(
    "hadamard"       => hadamard,
    "not"            => not,
    "cnot"           => cnot,
    "ccnot"          => ccnot,
    "phaseshift"     => phaseshift,
    "paulix"         => paulix,
    "pauliy"         => pauliy,
    "pauliz"         => pauliz,
    "swap"           => swap,
    "cswap"          => cswap,
    "sqrtswap"       => sqrtswap,
    "probe"          => probe,
    "measure"        => measure,
    "partialmeasure" => partialmeasure
)

const SYMBOLS = [("pi", pi), ("i", im)]

type Gate
    output::Int
    fn::Function
    inputs::Vector{Int}
    arguments::Vector{Any}
end

type InitialState
    output::Int
    state::QuantumState
end

function skipwhitespace(stream)
    pos = 1
    while pos <= length(stream) && isspace(stream[pos])
        pos += 1
    end
    return stream[pos:end]
end

function gateoutput(stream)
    stream = skipwhitespace(stream)

    if isempty(stream) || !isdigit(stream[1])
        return stream, Nullable{Int}()
    end

    lastindex = findfirst(stream, ':')
    return stream[lastindex + 1:end], tryparse(Int, stream[1:lastindex - 1])
end

function gatename(stream)
    stream = skipwhitespace(stream)

    if isempty(stream) || !isalpha(stream[1])
        return stream, Nullable{ASCIIString}()
    end

    lastindex = findfirst(c -> c == '[' || c == '(', stream)

    if lastindex == 0
        return stream, Nullable{ASCIIString}()
    end

    return stream[lastindex:end], Nullable(stream[1:lastindex - 1])
end

function inputs(stream)
    stream = skipwhitespace(stream)

    if isempty(stream) || stream[1] != '['
        return stream, Nullable{Vector{Int}}()
    end

    lastindex = findfirst(stream, ']')
    splitinputs = split(stream[2:lastindex - 1], ',', keep = false)

    if isempty(splitinputs)
        return stream[lastindex + 1:end], Nullable{Vector{Int}}()
    end

    inputs = map(x -> parse(Int, strip(x)), splitinputs)::Vector{Int}
    return stream[lastindex + 1:end], Nullable(inputs)
end

function tryparse_floatint(arg::AbstractString)
    tryint = tryparse(Int, arg)
    return !isnull(tryint) ? tryint : tryparse(Float64, arg)
end

function tryparse_symbol(special)
    symbol = SYMBOLS[1][1]
    symbolval = SYMBOLS[1][2]
    hassymbol = contains(special, symbol)

    if !hassymbol
        for i = 2:length(SYMBOLS)
            if contains(special, SYMBOLS[i][1])
                symbol = SYMBOLS[i][1]
                symbolval = SYMBOLS[i][2]
                hassymbol = true
                break
            end
        end
    end

    if !hassymbol
        return Nullable{Number}()
    end

    multiplier = split(special, symbol, keep = false)
    if isempty(multiplier)
        return Nullable(symbolval)
    end

    if length(multiplier) != 1
        error("Syntactic error in an argument expression")
    end

    trynum = tryparse_floatint(multiplier[1])
    return isnull(trynum) ? Nullable{Number}() : Nullable(get(trynum) * symbolval)
end

# A very primitive expression parser. It supports *one of* +-*/,
# the special constants i and pi, floating point and integer numbers
function tryparse_expression(expression)
    opindex = max(findfirst(expression, '+'),
                  findfirst(expression, '-'),
                  findfirst(expression, '*'),
                  findfirst(expression, '/'))

    if opindex == 0
        return Nullable{Number}()
    end

    # Now we only support a single operator
    op = expression[opindex]
    splitargs = split(expression, op)
    if length(splitargs) != 2
        return Nullable{Number}()
    end

    sides = []
    for arg in splitargs
        arg = strip(arg)
        result = tryparse_floatint(arg)

        if isnull(result)
            result = tryparse_symbol(arg)
        end

        if isnull(result)
            return Nullable{Number}()
        end

        push!(sides, get(result))
    end

    if op == '+'
        return Nullable(sides[1] + sides[2])
    elseif op == '-'
        return Nullable(sides[1] - sides[2])
    elseif op == '*'
        return Nullable(sides[1] * sides[2])
    elseif op == '/'
        return Nullable(sides[1] / sides[2])
    end

    return Nullable{Number}()
end

function argument(arg)
    arg = strip(arg)

    trynum = tryparse_floatint(arg)
    if !isnull(trynum)
        return get(trynum)
    end

    # Now try and parse it as a supported expression
    tryexpr = tryparse_expression(arg)

    # If nothing worked, we'll treat it as a string
    return isnull(tryexpr) ? arg : get(tryexpr)
end

function arguments(stream)
    stream = skipwhitespace(stream)

    if isempty(stream) || stream[1] != '('
        return stream, Nullable{Vector{Any}}()
    end

    lastindex = findfirst(stream, ')')
    splitargs = split(stream[2:lastindex - 1], ',', keep = false)

    if isempty(splitargs)
        return stream[lastindex + 1:end], Nullable{Vector{Any}}()
    end

    arguments = convert(Vector{Any}, map(argument, splitargs))::Vector{Any}
    return stream[lastindex + 1:end], Nullable(arguments)
end

function qubit(qs::AbstractString)
    if lowercase(qs) == "bell"
        return BELL_STATE
    end
    error("Unsupported qubit state")
end

function qubit(qi::Int)
    if qi == 0
        return QUBIT0
    elseif qi == 1
        return QUBIT1
    end
    error("Unsupported qubit state")
end

function parsenetlist(filename)
    file = open(filename)
    netlist = strip(readall(file))::ASCIIString
    close(file)

    initial::Vector{InitialState} = []
    gates::Vector{Gate} = []

    while !isempty(netlist)
        netlist, output = gateoutput(netlist)
        if isnull(output)
            error("No output index found for gate in the netlist")
        end

        netlist, name = gatename(netlist)::Tuple{ASCIIString, Nullable{ASCIIString}}
        if isnull(name)
            error("No name found for gate in the netlist")
        end

        netlist, inputs = Netlist.inputs(netlist)
        netlist, arguments = Netlist.arguments(netlist)

        if get(name) == "qstate"
            if isnull(arguments)
                error("No qubits given for an initial quantum state")
            end
            push!(initial, InitialState(get(output), fromstates(map(qubit, get(arguments))...)))
        else
            if !haskey(FUNCTION_MAP, get(name))
                error("Unsupported quantum gate detected: ", get(name))
            end
            fn = FUNCTION_MAP[get(name)]

            actualargs = isnull(arguments) ? [] : get(arguments)
            push!(gates, Gate(get(output), fn, get(inputs), actualargs))
        end
    end
    outputs = fill(Nullable{Any}(), length(initial) + length(gates))
    return (initial, gates, outputs)
end

end
