module QSpice

include("bitops.jl")
include("state.jl")
include("gates.jl")
include("parser.jl")
include("netlist.jl")

using QSpice.BitOps, QSpice.Gates, QSpice.Netlist, QSpice.State, QSpice.Parser

export simulate

function simulatefile(filename)
    f = open(filename)
    s = readstring(f)
    close(f)

    circuit = readnetlist(s)
    simulate(circuit)
end

end
