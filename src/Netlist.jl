module Netlist

using Iterators

using QSpice.State, QSpice.Gates

export parse, read_gates

function read_gates(filename::AbstractString)
    file = open(filename);
    gate_strings = imap(l -> split(l), eachline(file)) |> collect
    close(file)
    return vcat(gate_strings...)
end

end
