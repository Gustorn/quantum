include("QSpice.jl")

using QSpice, QSpice.Gates, QSpice.State, QSpice.Util, QSpice.Netlist
using Benchmarks


simulate("data/test.qnl")
simulate("data/test.qnl")

println("Manual reproduction of netlist")
q = from_states(QUBIT_0, QUBIT_1, QUBIT_1, QUBIT_0)
q = hadamard(q, 1)
q = pauli_y(q, 1)
q = swap(q, 1, 2)

probe(q, "Probe after swap")

q = hadamard(q, 3)
q = cnot(q, 1, 4,)
probe(q, "Probe after CNOT")

q = State.from_states(QUBIT_1, QUBIT_0, BELL_STATE)
isqrt = 1.0 / sqrt(2)

q = hadamard(q, 2)
q = cnot(q, 3, 2)
q = pauli_y(q, 2)
q = hadamard(q, 1)

measurement, posterior = partial_measure(q, 2)

print_bases(q)
@show measurement
print_bases(posterior)

println()
println()

println("Original state:")
print_bases(q)
println("Measuring bits 2 and 3")
m2, p2 = partial_measure(q, 2, 3)
@show m2
print_bases(p2)


#=
q = from_states(BELL_STATE, QUBIT_0, QUBIT_1, QUBIT_1, BELL_STATE, QUBIT_1, QUBIT_0)
@benchmark cnot(q, rand(1:q.bits), rand(1:q.bits))


#q1 = from_states(QUBIT_0, QUBIT_1, QUBIT_1)
#
#@benchmark swap(q1, rand(1:q1.bits), rand(1:q1.bits))

=#
