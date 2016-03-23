module Parser

using Iterators
using QSpice.Gates, QSpice.State, QSpice.Netlist

export netlist, flush

type GateStub
    behavior::Function
    arguments::Vector{Any}
    qedge::Vector{Int}
    bedge::Vector{Int}
end

const FUNCTION_MAP = Dict{ASCIIString, Function}(
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
    "superposition"  => superposition
)

function skipspace(s)
    pos = 1
    while !isempty(s[pos:end]) && isspace(s[pos])
        pos += 1
    end
    return s[pos:end]
end

function consume(s, token)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    return s
end

function expect(s, token, message)
    if startswith(s, token)
        return skipspace(s[length(token) + 1:end])
    end
    error(message)
end

function number(s)
    tryint = tryparse(Int ,s)
    return isnull(tryint) ? tryparse(Float64, s) : tryint
end

function namedfunction(s)
    for (name, fn) in FUNCTION_MAP
        if startswith(s, name)
            s = s[length(name) + 1:end]
            return skipspace(s), Nullable{Function}(fn)
        end
    end
    return s, Nullable{Function}()
end

function index(s)
    last = findfirst(s, ':')

    if last == 0 || isempty(s) || !isdigit(s[1])
        return s, Nullable{Int}()
    end

    tryint = tryparse(Int, s[1:last - 1])
    if isnull(tryint)
        return s, Nullable{Int}()
    end

    s = expect(s[last:end], ':', "Gate separator not found in netlist")
    return skipspace(s), Nullable(get(tryint))
end

function taggedindex(s, tag)
    if startswith(s, tag)
        last = findfirst(x -> x == ',' || x == ')', s)
        return skipspace(s[last:end]),
               tryparse(Int, s[length(tag) + 1:last - 1])
    end
    return s, Nullable{Int}()
end

function fallbackargs(s)
    last = findfirst(x -> x == ',' || x == ')', s)
    trynum = number(strip(s[1:last - 1]))
    return skipspace(s[last:end]), isnull(trynum) ? Nullable(s[1:last - 1]) : trynum
end

function simplefn(s)
    s, fn = namedfunction(s)

    if isnull(fn)
        return s, Nullable{Tuple{Function, Vector{Any}}}()
    end

    s = expect(s, '(', "Every function must have at least one argument")

    arguments = Vector{Any}([])
    while !startswith(s, ')')
        s, fnarg = simplefn(s)
        if !isnull(fnarg)
            fnchain = Vector{Tuple{Function, Vector{Any}}}()
            push!(fnchain, get(fnarg))

            while startswith(s, "|>")
                s, chained = simplefn(skipspace(s[3:end]))
                if isnull(chained)
                    error("There was a composition operator, but the following function was invalid")
                end
                push!(fnchain, get(chained))
            end
            push!(arguments, fnchain)
            s = consume(s, ',')
        else
            s, fba = fallbackargs(s)
            s = consume(s, ',')
            push!(arguments, get(fba))
        end
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments))
end

function gatefn(s)
    s, fn = namedfunction(s)

    if isnull(fn)
        return s, Nullable{Tuple{Function, Vector{Any}, Vector{Int}, Vector{Int}}}()
    end

    s = expect(s, '(', "Every function must have at least one argument")

    arguments = Vector{Any}([])
    qedges = Vector{Int}([])
    bedges = Vector{Int}([])

    while !startswith(s, ')')
        s, fnarg = simplefn(s)
        if !isnull(fnarg)
            fnchain = Vector{Tuple{Function, Vector{Any}}}()
            push!(fnchain, get(fnarg))

            while startswith(s, "|>")
                s, chained = simplefn(skipspace(s[3:end]))
                if isnull(chained)
                    error("There was a composition operator, but the following function was invalid")
                end
                push!(fnchain, get(chained))
            end
            push!(arguments, fnchain)
            s = consume(s, ',')
            continue
        end

        s, qedge = taggedindex(s, 'Q')
        if !isnull(qedge)
            s = consume(s, ',')
            push!(qedges, get(qedge))
            continue
        end

        s, bedge = taggedindex(s, 'B')
        if !isnull(bedge)
            s= consume(s, ',')
            push!(bedges, get(bedge))
            continue
        end

        s, fba = fallbackargs(s)
        s = consume(s, ',')
        push!(arguments, get(fba))
    end
    return skipspace(s[2:end]), Nullable((get(fn), arguments, qedges, bedges))
end

function qubit(q::AbstractString)
    if q == "0"
        return QUBIT0
    elseif q == "1"
        return QUBIT1
    elseif lowercase(q) == "bell"
        return BELL_STATE
    elseif lowercase(q) == "rand"
        s = [rand() + rand() * im, rand() + rand() * im]
        n = sqrt(sumabs2(s))
        s = s .* (1.0 / n)
        qs = QuantumState(s, 1)
        return qs
    end
    error("Unsupported qubit found in netlist")
end

function bit(b::AbstractString)
    if b == "0"
        return 0
    elseif b == "1"
        return 1
    end
    error("Classical bits can only take a value of 0 or 1")
end

function gate(s)
    rollback = s
    s, idx = index(s)
    if isnull(idx)
        return rollback, Nullable{Tuple{Int, GateStub}}()
    end

    s, tryfn = gatefn(s)
    if isnull(tryfn)
        return rollback, Nullable{Tuple{Int, GateStub}}()
    end

    fn, args, qedge, bedge = get(tryfn)
    stub = GateStub(fn, args, qedge, bedge)
    return skipspace(s), Nullable((get(idx), stub))
end

function qstate(s)
    rollback = s
    s, idx = index(s)
    if isnull(idx)
        return rollback, Nullable{Tuple{Int, QuantumState}}()
    end

    nameend = findfirst(s, '(')
    argend = findfirst(s, ')')

    name = s[1:nameend - 1]
    argstring = s[nameend + 1:argend - 1]

    if name == "qstate"
        args = map(x -> qubit(strip(x)), split(argstring, ',', keep = false))
        return skipspace(s[argend + 1:end]), Nullable((get(idx), fromstates(args...)))
    end
    return rollback, Nullable{Tuple{Int, QuantumState}}()
end

function bstate(s)
    rollback = s
    s, idx = index(s)
    if isnull(idx)
        return rollback, Nullable{Tuple{Int, Vector{Int}}}()
    end

    nameend = findfirst(s, '(')
    argend = findfirst(s, ')')

    name = s[1:nameend - 1]
    argstring = s[nameend + 1:argend - 1]

    if name == "bstate"
        args = map(x -> bit(strip(x)), split(argstring, ',', keep = false))
        return skipspace(s[argend + 1:end]), Nullable((get(idx), args))
    elseif name == "randbstate"
        nbits = parse(Int, strip(argstring))
        return skipspace(s[argend + 1:end]), Nullable((get(idx), [rand(0:1) for _ in 1:nbits]))
    end
    return rollback, Nullable{Tuple{Int, Vector{Int}}}()
end

function rawelements(s)
    stubs::Vector{Tuple{Int, GateStub}} = []
    qstates::Vector{Tuple{Int, QuantumState}} = []
    bstates::Vector{Tuple{Int, Vector{Int}}} = []
    while !isempty(s)
        s, trygate = gate(s)
        if !isnull(trygate)
            push!(stubs, get(trygate))
            continue
        end

        s, tryqstate = qstate(s)
        if !isnull(tryqstate)
            push!(qstates, get(tryqstate))
            continue
        end

        s, trybstate = bstate(s)
        if !isnull(trybstate)
            push!(bstates, get(trybstate))
            continue
        end
        error("Invalid netlist")
    end
    return stubs, qstates, bstates
end

function addtodag(key::Int, value::GateStub, dag::Dict)
    if !haskey(dag, key)
        dag[key] = GateNode(value.behavior, value.arguments, Nullable(), Nullable(), [], [])
    end
end

function addtodag(key::Int, value::QuantumState, dag::Dict)
    if !haskey(dag, key)
        dag[key] = QuantumStateNode(value, [])
    end
end

function addtodag(key::Int, value::Vector{Int}, dag::Dict)
    if !haskey(dag, key)
        dag[key] = BitStateNode(value, [])
    end
end

function netlist(s)
    gatestubs, qstubs, bstubs = rawelements(s)
    merged = [i[1] => i[2] for i in chain(gatestubs, qstubs, bstubs)]
    dag = Dict()


    for (key, value) in gatestubs
        addtodag(key, value, dag)

        for qi in value.qedge
            addtodag(qi, merged[qi], dag)
            push!(dag[qi].qedge, dag[key])
        end

        for bi in value.bedge
            addtodag(bi, merged[bi], dag)
            push!(dag[bi].bedge, dag[key])
        end
    end

    return map(x -> x[2], dag)
end

function flush(dag)
    for node in dag
        if isa(node, QuantumStateNode)
            Netlist.broadcast(node)
        elseif isa(node, BitStateNode)
            Netlist.broadcast(node)
        end
    end
end
end
