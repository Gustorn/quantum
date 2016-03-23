module Netlist

using Match

using QSpice.State, QSpice.Gates

export GateNode, QuantumStateNode, BitStateNode, broadcast, accept

abstract GateType

type Nullary <: GateType end

type Unary{T} <: GateType
    param::T
end

type Binary{T, U} <: GateType
    param1::T
    param2::U
end

type Ternary{T, U, V} <: GateType
    param1:T
    param2:U
    param3:V
end

type Dynamic{T} <: GateType
    params::Vector{T}
end

type GateNode{T <: GateType, U <: GateType, V <: GateType}
    behavior::Function

    arguments::T
    quantuminputs::U
    bitinputs::V

    qedge::Vector{GateNode}
    bedge::Vector{GateNode}
end

type QuantumStateNode
    state::QuantumState
    qedge::Vector{GateNode}
end

type BitStateNode
    bits::Vector{Int}
    bedge::Vector{GateNode}
end

const QUANTUM_ONLY = [qidentity, hadamard, not, cnot,
                      ccnot, swap, cswap, sqrtswap,
                      paulix, pauliy, pauliz, phaseshift,
                      measure, partialmeasure, probe, superposition]

const CLASSICAL_ONLY = [probe]

const HYBRID = [choose1]

function gatenode(behavior::Function, args::Vector{Any})
    nullquantum() = Unary{QuantumState}(Nullable{QuantumState}())
    nullbits() = Unary{Vector{Int}}(Nullable{Vector{Int}}())

    @match behavior begin
        qidentity      => GateNode(behavior, Nullary(),                       nullquantum(), Nullary(),  [], [])
        hadamard       => GateNode(behavior, Unary{Int}(args...),             nullquantum(), Nullary(),  [], [])
        not            => GateNode(behavior, Unary{Int}(args...),             nullquantum(), Nullary(),  [], [])
        cnot           => GateNode(behavior, Binary{Int, Int}(args...),       nullquantum(), Nullary(),  [], [])
        ccnot          => GateNode(behavior, Ternary{Int, Int, Int}(args...), nullquantum(), Nullary(),  [], [])
        swap           => GateNode(behavior, Binary{Int, Int}(args...),       nullquantum(), Nullary(),  [], [])
        cswap          => GateNode(behavior, Ternary{Int, Int, Int}(args...), nullquantum(), Nullary(),  [], [])
        sqrtswap       => GateNode(behavior, Binary{Int, Int}(args...),       nullquantum(), Nullary(),  [], [])
        paulix         => GateNode(behavior, Unary{Int}(args...),             nullquantum(), Nullary(),  [], [])
        pauliy         => GateNode(behavior, Unary{Int}(args...),             nullquantum(), Nullary(),  [], [])
        pauliz         => GateNode(behavior, Unary{Int}(args...),             nullquantum(), Nullary(),  [], [])
        phaseshift     => GateNode(behavior, Binary{Int, Float64}(args...),   nullquantum(), Nullary(),  [], [])
        measure        => GateNode(behavior, Nullary(),                       nullquantum(), Nullary(),  [], [])
        partialmeasure => GateNode(behavior, Dynamic{Int},                    nullquantum(), Nullary(),  [], [])
        probe          => GateNode(behavior, Unary{AbstractString}(args...),  nullquantum(), nullbits(), [], [])
        superposition  => GateNode(behavior, Nullary(),                       Dynamic{QuantumState}(Nullable()), Nullary(), [], [])
    end
end

function broadcast(gate::GateNode, qstate::QuantumState)
    for q in gate.qedge
        accept(q, qstate)
    end
end

function broadcast(gate::QuantumStateNode)
    for q in gate.qedge
        accept(q, gate.state)
    end
end

function broadcast(gate::GateNode, bstate::Vector{Int})
    for b in gate.bedge
        accept(b, bstate)
    end
end

function broadcast(gate::BitStateNode)
    for b in gate.bedge
        accept(b, gate.bits)
    end
end

function broadcast(gate::GateNode, output::Tuple{QuantumState, Vector{Int}})
    broadcast(gate, output[1])
    broadcast(gate, output[2])
end

function broadcast(gate::GateNode, x::Void) end

function executeready(gate::GateNode)
    if gate.behavior in HYBRID && !isnull(gate.qinput) && !isnull(gate.binput)
        broadcast(gate, gate.behavior(get(gate.qinput), get(gate.binput), gate.arguments...))
        gate.qinput = Nullable()
        gate.binput = Nullable()
    end
end

function accept(gate::GateNode, qstate::QuantumState, states::QuantumState...)
    if gate.behavior in QUANTUM_ONLY
        broadcast(gate, gate.behavior(qstate, states..., gate.arguments...))
    end
end

function accept(gate::GateNode, qstate::QuantumState)
    if gate.behavior in QUANTUM_ONLY
        broadcast(gate, gate.behavior(qstate, gate.arguments...))
        gate.binput = Nullable()
    else
        gate.qinput = Nullable(qstate)
        executeready(gate)
    end
end

function accept(gate::GateNode, bstate::Vector{Int})
    if gate.behavior in CLASSICAL_ONLY
        broadcast(gate, gate.behavior(bstate, gate.arguments...))
        gate.qinput = Nullable()
    else
        gate.binput = Nullable(bstate)
        executeready(gate)
    end
end

function accept(gate::GateNode, qstate::QuantumState, bstate::Vector{Int})
    if gate.behavior in HYBRID
        broadcast(gate, gate.behavior(get(gate.qinput), get(gate.binput), gate.arguments...))
        gate.qinput = Nullable()
        gate.binput = Nullable()
    else
        accept(gate, qstate)
        accept(gate, bstate)
    end
end
end
