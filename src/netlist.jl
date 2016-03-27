module Netlist

using QSpice.State, QSpice.Gates, QSpice.Parser

export QuantumCircuit, readnetlist, simulate

# Any subtype of this may be a node in the DAG representing the
# quantum circuit
abstract QuantumCircuitNode

type GateNode <: QuantumCircuitNode
    fn::Function
    args::Vector{Any}

    qinputs::Vector{QuantumState}
    binputs::Vector{Vector{Int}}

    quantumedges::Vector{Int}
    bitedges::Vector{Int}
end

type QuantumNode <: QuantumCircuitNode
    state::QuantumState
    edges::Vector{Int}
end

type BitNode <: QuantumCircuitNode
    bits::Vector{Int}
    edges::Vector{Int}
end

type QuantumCircuit
    gatenodes::Dict{Int, GateNode}
    qnodes::Dict{Int, QuantumNode}
    bnodes::Dict{Int, BitNode}
end

# The list of nodes that don't need to wait for additional inputs
# in case of a quantum state change
const QUANTUM_PASSTHROUGH = [Gates.identity, hadamard, not, cnot,
                             ccnot, swap, cswap, sqrtswap, qft,
                             paulix, pauliy, pauliz, phaseshift,
                             measure, partialmeasure, probe]

# The list of nodes that don't need to wait for additional inputs
# in case of classical bit state change
const BIT_PASSTHROUGH = [probe]

function readnetlist(s)
    function addnotfound(gatenodes, index, gate)
        if !haskey(gatenodes, index)
            gatenodes[index] = GateNode(gate.fn, gate.args, [], [], [], [])
        end
    end

    gates, qstates, bstates = parsenetlist(s)

    gatenodes = Dict{Int, GateNode}()

    # First initialize quantum and classical states
    qnodes = [i[1] => QuantumNode(i[2], []) for i in qstates]
    bnodes = [i[1] => BitNode(i[2], []) for i in bstates]

    # Construct the DAG by iterating the gates and adding all necessary
    # edges (representing a dependency on either a quantum or classical state change)
    # to all of the nodes
    for (index, gate) in gates
        addnotfound(gatenodes, index, gate)

        for qi in gate.quantuminputs
            if haskey(qnodes, qi)
                push!(qnodes[qi].edges, index)
                continue
            end
            addnotfound(gatenodes, qi, gates[qi])
            push!(gatenodes[qi].quantumedges, index)
        end

        for bi in gate.bitinputs
            if haskey(bnodes, bi)
                push!(bnodes[bi].edges, index)
                continue
            end
            addnotfound(gatenodes, bi, gates[bi])
            push!(gatenodes[bi].bitedges, index)
        end
    end

    return QuantumCircuit(gatenodes, qnodes, bnodes)
end

# Broadcasts a quantum state change to all nodes connected by a quantum edge
function broadcast(netlist::QuantumCircuit, gatenode::GateNode, qstate::QuantumState)
    for qedge in gatenode.quantumedges
        accept(netlist, netlist.gatenodes[qedge], qstate)
    end
end

# Broadcasts a classical bit state change to all nodes connected by a bit edge
function broadcast(netlist::QuantumCircuit, gatenode::GateNode, bstate::Vector{Int})
    for bedge in gatenode.bitedges
        accept(netlist, netlist.gatenodes[bedge], bstate)
    end
end

# Broadcasts both a quantum and classical state change
function broadcast(netlist::QuantumCircuit, gatenode::GateNode,
                   states::Tuple{QuantumState, Vector{Int}})
    broadcast(netlist, gatenode, states[1])
    broadcast(netlist, gatenode, states[2])
end

function broadcast(netlist::QuantumCircuit, gatenode::GateNode, v::Void)
    # Ordinarily only probes should call this function
end

# Used by gates that have multiple inputs. It checks if the gate is ready to execute,
# and if it is, performs the operation.
# TODO(gustorn): this only supports unordered arguments right now, which is a pretty big limiation
# (and a wonder that it even works in case of superposition). Maybe make edges between inputs/outputs
# and not gates?
function executeready(netlist::QuantumCircuit, gatenode::GateNode)
    if gatenode.fn == superposition && length(gatenode.qinputs) == 2
        broadcast(netlist, gatenode, gatenode.fn(gatenode.qinputs..., gatenode.args...))
        empty!(gatenode.qinputs)
    elseif gatenode.fn == choose1 && !isempty(gatenode.qinputs) && !isempty(gatenode.binputs)
        broadcast(netlist, gatenode, gatenode.fn(gatenode.qinputs...,
                                                 gatenode.binputs...,
                                                 convert(Vector{Vector{Tuple{Function, Vector{Any}}}}, gatenode.args)))
        empty!(gatenode.qinputs)
        empty!(gatenode.binputs)
    end
end

function accept(netlist::QuantumCircuit, gatenode::GateNode, qstate::QuantumState)
    if gatenode.fn in QUANTUM_PASSTHROUGH
        broadcast(netlist, gatenode, gatenode.fn(qstate, gatenode.args...))
        return
    end
    push!(gatenode.qinputs, qstate)
    executeready(netlist, gatenode)
end

function accept(netlist::QuantumCircuit, gatenode::GateNode, bstate::Vector{Int})
    if gatenode.fn in BIT_PASSTHROUGH
        broadcast(netlist, gatenode, gatenode.fn(bstate, gatenode.args...))
        return
    end
    push!(gatenode.binputs, bstate)
    executeready(netlist, gatenode)
end

# Simulates the netlist by flushing all initial qubits and classical bitstates.
# If the DAG is correct, this will execute every reachable gate in the circuit
function simulate(netlist::QuantumCircuit)
    for qn in values(netlist.qnodes)
        for edge in qn.edges
            accept(netlist, netlist.gatenodes[edge], qn.state)
        end
    end
    for bn in values(netlist.bnodes)
        for edge in bn.edges
            accept(netlist, netlist.gatenodes[edge], bn.bits)
        end
    end
end

end
