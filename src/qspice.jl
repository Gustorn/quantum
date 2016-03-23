module QSpice

using Iterators.imap

include("bitops.jl")
include("state.jl")
include("gates.jl")
include("netlist.jl")
include("parser.jl")

using QSpice.BitOps, QSpice.Gates, QSpice.Netlist, QSpice.State, QSpice.Parser

export simulate

function simulate(filename)
    f = open(filename)
    s = readall(f)
    close(f)

    circuit = netlist(s)
    Parser.flush(circuit)
end

end
