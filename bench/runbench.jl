include("../src/QSpice.jl")

using Benchmarks

using QSpice.Gates, QSpice.State

s = from_states(BELL_STATE, QUBIT_0, QUBIT_1,
                QUBIT_0, QUBIT_1, QUBIT_1, BELL_STATE,
                BELL_STATE, BELL_STATE)

cn = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

Benchmarks.@benchmarkable(hadamard!,         nothing, hadamard(s, rand(1:s.bits)),                              nothing)
Benchmarks.@benchmarkable(not!,              nothing, not(s, rand(1:s.bits)),                                   nothing)
Benchmarks.@benchmarkable(cnot!,             nothing, cnot(s, rand(1:s.bits), rand(1:s.bits)),                  nothing)
Benchmarks.@benchmarkable(ccnot!,            nothing, ccnot(s, rand(1:s.bits), rand(1:s.bits), rand(1:s.bits)), nothing)
Benchmarks.@benchmarkable(swap!,             nothing, swap(s, rand(1:s.bits), rand(1:s.bits)),                  nothing)
Benchmarks.@benchmarkable(cswap!,            nothing, cswap(s, rand(1:s.bits), rand(1:s.bits), rand(1:s.bits)), nothing)
Benchmarks.@benchmarkable(sqrt_swap!,        nothing, sqrt_swap(s, rand(1:s.bits), rand(1:s.bits)),             nothing)
Benchmarks.@benchmarkable(phase_shift!,      nothing, phase_shift(s, rand(1:s.bits), rand(1:2pi)),              nothing)
Benchmarks.@benchmarkable(pauli_x!,          nothing, pauli_x(s, rand(1:s.bits)),                               nothing)
Benchmarks.@benchmarkable(pauli_y!,          nothing, pauli_y(s, rand(1:s.bits)),                               nothing)
Benchmarks.@benchmarkable(pauli_z!,          nothing, pauli_z(s, rand(1:s.bits)),                               nothing)
Benchmarks.@benchmarkable(measure!,          nothing, measure(s),                                               nothing)
Benchmarks.@benchmarkable(partial_measure!,  nothing, partial_measure(s, rand(1:s.bits), rand(1:s.bits),
                                                                         rand(1:s.bits), rand(1:s.bits),
                                                                         rand(1:s.bits), rand(1:s.bits)),       nothing)
Benchmarks.@benchmarkable(unitary_cnot!,     nothing, unitary(s, cn, 4, 3),                                     nothing)


benchmarks = [hadamard!, not!, cnot!, ccnot!, swap!, cswap!, sqrt_swap!, phase_shift!, pauli_x!,
              pauli_y!, pauli_z!, measure!, partial_measure!, unitary_cnot!]

for b in benchmarks
    println("Benchmark: ", b)
    result = Benchmarks.execute(b)
    println(result)
end
