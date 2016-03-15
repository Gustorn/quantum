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

function parse_one(gate_string::AbstractString)

end

function parse(filename::AbstractString)
    gate_strings = read_gates(filename)
    outputs = fill(null(), length(gate_strings))

end

end
