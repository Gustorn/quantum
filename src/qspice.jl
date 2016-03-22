module QSpice

using Iterators.imap

include("bitops.jl")
include("state.jl")
include("gates.jl")
include("netlist.jl")

using QSpice.BitOps, QSpice.Gates, QSpice.Netlist, QSpice.State

export simulate

# TODO(gustorn): handle the outputs of measurement gates correctly
function flush(gates::Vector{Gate}, outputs::Vector{Nullable{Any}})
    for g in gates
        inputs = map(x -> outputs[x], g.inputs)
        if !any(isnull, inputs) && isnull(outputs[g.output])
            outputs[g.output] = g.fn(imap(get, inputs)..., g.arguments...)
        end
    end
end

function simulate(filename)
    initial, gates, outputs = parsenetlist(filename)
    for s in initial
        outputs[s.output] = s.state
    end

    todo = count(isnull, outputs)
    prev = -1

    iteration = 1
    while todo != prev
        println("================\nIteration: ",
                iteration,
                "\n================")
        flush(gates, outputs)
        prev = todo
        todo = count(isnull, outputs)
        iteration += 1
    end
end

end
