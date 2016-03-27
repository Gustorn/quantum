# Important note when extending the parser: make sure to call skipspace on the returned
# stream whenever the parsing method is successful. Since netlist format doesn't care about
# whitespace, this simplifies every other algorithm using these building blocks.
module Parser

using QSpice.Gates, QSpice.State

export ParsedGate, parsenetlist

type ParsedGate
    fn::Function
    args::Vector{Any}
    quantuminputs::Vector{Int}
    bitinputs::Vector{Int}
end

parsedgate(t::Tuple{Function, Vector{Any}, Vector{Int}, Vector{Int}}) = ParsedGate(t[1], t[2], t[3], t[4])

const FUNCTION_MAP = Dict{ASCIIString, Function}(
    "superposition"  => superposition,
    "identity"       => Gates.identity,
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
    "partialmeasure" => partialmeasure,
    "choose1"        => choose1,
    "qft"            => qft,
)

function skipspace(s)
    pos = 1
    while !isempty(s[pos:end]) && isspace(s[pos])
        pos += 1
    end
    return s[pos:end]
end

# Skips a token if it was present in the input stream, otherwise
# does nothing
function skip(s, token)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    return s
end

# Consumes the stream up to and including the terminator value,
# parsing and returning the value consumed with the type passed
# as the second argument.
function consume(s, Int, terminator)
    last = findfirst(s, terminator)
    if last == 0
        return s, Nullable{Int}()
    end
    return skipspace(s[last + 1:end]), tryint(s[1:last - 1])
end

# Skips a token if it was present in the inputs stream, otherwise
# throws an error
function expect(s, token)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    error("Expected token ", token, " not found in netlist")
end

tryint(s) = tryparse(Int, s)
tryfloat(s) = tryparse(Float64, s)

# Parser a tagged index value, with the format of <tag><int>,
# where <tag> may be an arbitrary string and <int> an 64-bit
# integer value. The two pieces may NOT be separated by whitespace.
function taggedindex(s, tag)
    if !startswith(s, tag)
        return s, Nullable{Int}()
    end

    last = findfirst(x -> x == ',' || x == ')', s)
    if last == 0
        return s, Nullable{Int}()
    end

    i = tryint(s[length(tag) + 1:last - 1])
    if isnull(i)
        return s, Nullable{Int}()
    end

    return skipspace(s[last:end]), Nullable(get(i))
end

# Used as a fallback for argument parsing. Tries to parse the argument
# as an integer, then as a floating point number, and if both fail it returns
# the argument as a string.
function fallback(s)
    last = findfirst(x -> x == ',' || x == ')', s)
    if last == 0
        error("Fallback parser called outside of an argument list")
    end

    i = tryint(s[1:last - 1])
    if !isnull(i)
        return skipspace(s[last:end]), get(i)
    end

    f = tryfloat(s[1:last - 1])
    if !isnull(f)
        return skipspace(s[last:end]), get(f)
    end

    return skipspace(s[last:end]), s[1:last - 1]
end

function gatename(s)
    for (name, fn) in FUNCTION_MAP
        if startswith(s, name)
            return skipspace(s[length(name) + 1:end]), Nullable(fn)
        end
    end
    return s, Nullable{Function}()
end

# Parser a netlist function, but without quantum/bit state inputs. These
# are used as arguments to top-level functions. Example:
# choose1(Q1, B1, choose1(not(1)), not(1))
# ^full           ^simple ^simple  ^simple
function simplefn(s)
    rollback = s

    s, fn = gatename(s)
    if isnull(fn)
        return rollback, Nullable{Tuple{Function, Vector{Any}}}()
    end

    if isempty(s) || s[1] != '('
        return rollback, Nullable{Tuple{Function, Vector{Any}}}()
    end
    s = skip(s, '(')

    arguments = Vector{Any}([])
    while !startswith(s, ')')
        # First see if the argument in a function
        s, argfn = simplefn(s)
        if !isnull(argfn)
            # If it is, parse the full chain
            fnchain = Vector{Tuple{Function, Vector{Any}}}([])
            push!(fnchain, get(argfn))

            while startswith(s, "|>")
                s, chainfn = simplefn(skipspace(s[3:end]))
                push!(fnchain, get(chainfn))
            end
            push!(arguments, fnchain)
            s = skip(s, ',')
            continue
        end

        # If not, use the fallback parser
        s, fb = fallback(s)
        push!(arguments, fb)
        s = skip(s, ',')
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments))
end

# Parses a top level function
function fullfn(s)
    rollback = s
    s, fn = gatename(s)
    if isnull(fn)
        return rollback, Nullable{Tuple{Function, Vector{Any}, Vector{Int}, Vector{Int}}}()
    end

    if isempty(s) || s[1] != '('
        return rollback, Nullable{Tuple{Function, Vector{Any}, Vector{Int}, Vector{Int}}}()
    end
    s = skip(s, '(')

    arguments = Vector{Any}([])
    qedges = Vector{Int}([])
    bedges = Vector{Int}([])

    while !startswith(s, ')')
        # First see if the argument in a function
        s, argfn = simplefn(s)
        if !isnull(argfn)
            # If it is, parse the full chain
            fnchain = Vector{Tuple{Function, Vector{Any}}}([])
            push!(fnchain, get(argfn))

            while startswith(s, "|>")
                s, chainfn = simplefn(skipspace(s[3:end]))
                push!(fnchain, get(chainfn))
            end
            push!(arguments, fnchain)
            s = skip(s, ',')
            continue
        end

        # If it's not, see if it's a quantum connection
        s, qi = taggedindex(s, 'Q')
        if !isnull(qi)
            push!(qedges, get(qi))
            s = skip(s, ',')
            continue
        end

        # Then if it's a classical connection
        s, bi = taggedindex(s, 'B')
        if !isnull(bi)
            push!(bedges, get(bi))
            s = skip(s, ',')
            continue
        end

        # And as a last resort use the fallback parser
        s, fb = fallback(s)
        push!(arguments, fb)
        s = skip(s, ',')
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments, qedges, bedges))
end

# Tries to parse a single gate with the syntax: <index>:<name>(args...)
function gate(s)
    rollback = s
    s, index = consume(s, Int, ':')
    if isnull(index)
        return rollback, Nullable{Tuple{Int, ParsedGate}}()
    end

    s, parsed = fullfn(s)
    if isnull(parsed)
        return rollback, Nullable{Tuple{Int, ParsedGate}}()
    end

    return skipspace(s), Nullable((get(index), parsedgate(get(parsed))))
end

# Tries to parse a single quantum state with the syntax: <index>:qstate(qubit),
# where qubit is one of: 0, 1, bell, rand
function qstate(s)
    function qubit(q)
        if q == "0"
            return QUBIT0
        elseif q == "1"
            return QUBIT1
        elseif q == "bell"
            return BELL_STATE
        elseif q == "rand"
            # Generate a random qubit
            state = [rand([-1, 1]) * rand() + rand([-1, 1]) * rand() * im for _ in 1:2]
            # Then normalize it
            n = sqrt(sumabs2(state))
            state = state .* (1.0 / n)
            return QuantumState(state, 1)
        end
        error("Unsupported qubit \"", q, "\" in netlist")
    end

    rollback = s
    s, index = consume(s, Int, ':')
    if isnull(index)
        return rollback, Nullable{Tuple{Int, QuantumState}}()
    end

    if !startswith(s, "qstate")
        return rollback, Nullable{Tuple{Int, QuantumState}}()
    end

    s = skip(s, "qstate")
    s = expect(s, '(')
    last = findfirst(s, ')')
    states = map(x -> qubit(strip(x)), split(s[1:last - 1], ',', keep = false))::Vector{QuantumState}
    return skipspace(s[last + 1:end]), Nullable((get(index), fromstates(states...)))
end

# Tries to parse a single classical bit state with syntax: <index>:bstate(bit...)
# where bit is one of the following: 0, 1, rand
function bstate(s)
    function bit(b)
        if b == "0"
            return 0
        elseif b == "1"
            return 1
        elseif b == "rand"
            return rand(0:1)
        end
        error("A classical bit can only be 0 or 1")
    end

    rollback = s
    s, index = consume(s, Int, ':')
    if isnull(index)
        return rollback, Nullable{Tuple{Int, Vector{Int}}}()
    end

    if !startswith(s, "bstate")
        return rollback, Nullable{Tuple{Int, Vector{Int}}}()
    end

    s = skip(s, "bstate")
    s = expect(s, '(')
    last = findfirst(s, ')')
    bits = map(x -> bit(strip(x)), split(s[1:last - 1], ',', keep = false))::Vector{Int}
    return skipspace(s[last + 1:end]), Nullable((get(index), bits))
end

# Tries to parse the whole netlist
function parsenetlist(s)
    gates = Dict{Int, ParsedGate}()
    qstates = Dict{Int, QuantumState}()
    bstates = Dict{Int, Vector{Int}}()
    while !isempty(s)
        s, trygate = gate(s)
        if !isnull(trygate)
            idx, g = get(trygate)
            gates[idx] = g
            continue
        end

        s, tryqstate = qstate(s)
        if !isnull(tryqstate)
            idx, q = get(tryqstate)
            qstates[idx] = q
            continue
        end

        s, trybstate = bstate(s)
        if !isnull(trybstate)
            idx, b = get(trybstate)
            bstates[idx] = b
            continue
        end
        error("Invalid syntax in netlist, starting at: ", s)
    end
    return gates, qstates, bstates
end

end
