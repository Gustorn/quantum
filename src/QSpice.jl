module QSpice

using Iterators.imap

include("Util.jl")
include("State.jl")
include("Gates.jl")
include("Netlist.jl")

using QSpice.State, QSpice.Gates, QSpice.Util, QSpice.Netlist

export simulate

function count_unfinished(outputs::Vector{Nullable{QuantumState}})
    return count(isnull, outputs)
end

function flush(gates::Vector{Gate}, outputs::Vector{Nullable{QuantumState}})
    for g in gates
        inputs = map(x -> outputs[x], g.inputs)
        if !any(isnull, inputs) && isnull(outputs[g.index])
            outputs[g.index] = g.fn(imap(get, inputs)..., g.arguments...)
        end
    end
end

function simulate{S<:AbstractString}(filename::S)
    initial_states, gates, outputs = parse_netlist(filename)
    for s in initial_states
        outputs[s.index] = s.state
    end

    prev_unfinished = count_unfinished(outputs)
    curr_unfinished = -1

    iteration = 1
    while prev_unfinished != curr_unfinished
        println("================\nIteration: ",
                iteration,
                "\n================")
        flush(gates, outputs)
        prev_unfinished = curr_unfinished
        curr_unfinished = count_unfinished(outputs)
        iteration += 1
    end
end

end
