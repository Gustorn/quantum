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
    "identity"       => qidentity,
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

function skip(s, token)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    return s
end

function consume(s, Int, terminator)
    last = findfirst(s, terminator)
    if last == 0
        return s, Nullable{Int}()
    end
    return skipspace(s[last + 1:end]), tryint(s[1:last - 1])
end

function expect(s, token)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    error("Expected token ", token, " not found in netlist")
end

tryint(s) = tryparse(Int, s)
tryfloat(s) = tryparse(Float64, s)

function taggedindex(s, tag)
    last = findfirst(x -> x == ',' || x == ')', s)
    if last == 0
        return s, Nullable{Int}()
    end
    if !startswith(s, tag)
        return s, Nullable{Int}()
    end

    i = tryint(s[length(tag) + 1:last - 1])
    if isnull(i)
        return s, Nullable{Int}()
    end

    return skipspace(s[last:end]), Nullable(get(i))
end

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
        s, argfn = fullfn(s)
        if !isnull(argfn)
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

        s, fb = fallback(s)
        push!(arguments, fb)
        s = skip(s, ',')
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments))
end

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
        s, argfn = simplefn(s)
        if !isnull(argfn)
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

        s, qi = taggedindex(s, 'Q')
        if !isnull(qi)
            push!(qedges, get(qi))
            s = skip(s, ',')
            continue
        end

        s, bi = taggedindex(s, 'B')
        if !isnull(bi)
            push!(bedges, get(bi))
            s = skip(s, ',')
            continue
        end

        s, fb = fallback(s)
        push!(arguments, fb)
        s = skip(s, ',')
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments, qedges, bedges))
end

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

function qstate(s)
    function qubit(q)
        if q == "0"
            return QUBIT0
        elseif q == "1"
            return QUBIT1
        elseif q == "bell"
            return BELL_STATE
        elseif q == "rand"
            state = [rand([-1, 1]) * rand() + rand([-1, 1]) * rand() * im for _ in 1:2]
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
