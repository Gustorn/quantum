include("QSpice.jl")

using QSpice.Gates, QSpice.State, QSpice.Util

q = State.from_states(BELL_STATE, QUBIT_1, QUBIT_1)
isqrt = 1.0 / sqrt(2)

q = hadamard(q, 2)
#q = cnot(q, 3, 2)
q = pauli_y(q, 2)
q = hadamard(q, 1)

#probe(q)
#println("Measurement:", measure(q))

measurement, posterior = partial_measure(q, 2)

print_bases(q)
@show measurement
print_bases(posterior)

println()
println()

println("Original state:")
print_bases(q)
println("Measuring bits 1 and 3")
m2, p2 = partial_measure(q, 2, 3)
@show m2
print_bases(p2)
