include("QSpice.jl")

using QSpice.Gates, QSpice.State, QSpice.Util

q = State.from_states(QUBIT_0, QUBIT_1, QUBIT_1)
isqrt = 1.0 / sqrt(2)

q = hadamard(q, 2)
q = pauli_y(q, 1)
q = hadamard(q, 1)

#probe(q)
#println("Measurement:", measure(q))
#probe(BELL_STATE)

measurement, posterior = partial_measure(q, 2)

print_bases(q)
@show measurement
print_bases(posterior)

