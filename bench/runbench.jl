include("../src/qspice.jl")

using Benchmarks

using QSpice.Gates, QSpice.State

s = fromstates(BELL_STATE, QUBIT0, QUBIT1,
               QUBIT0, QUBIT1, QUBIT1, BELL_STATE,
               BELL_STATE, BELL_STATE)

cn = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

Benchmarks.@benchmarkable(hadamard!,         nothing, hadamard(s, rand(1:s.bits)),                              nothing)
Benchmarks.@benchmarkable(not!,              nothing, not(s, rand(1:s.bits)),                                   nothing)
Benchmarks.@benchmarkable(cnot!,             nothing, cnot(s, rand(1:s.bits), rand(1:s.bits)),                  nothing)
Benchmarks.@benchmarkable(ccnot!,            nothing, ccnot(s, rand(1:s.bits), rand(1:s.bits), rand(1:s.bits)), nothing)
Benchmarks.@benchmarkable(swap!,             nothing, swap(s, rand(1:s.bits), rand(1:s.bits)),                  nothing)
Benchmarks.@benchmarkable(cswap!,            nothing, cswap(s, rand(1:s.bits), rand(1:s.bits), rand(1:s.bits)), nothing)
Benchmarks.@benchmarkable(sqrtswap!,         nothing, sqrtswap(s, rand(1:s.bits), rand(1:s.bits)),              nothing)
Benchmarks.@benchmarkable(phaseshift!,       nothing, phaseshift(s, rand(1:s.bits), rand(1:2pi)),               nothing)
Benchmarks.@benchmarkable(paulix!,           nothing, paulix(s, rand(1:s.bits)),                                nothing)
Benchmarks.@benchmarkable(pauliy!,           nothing, pauliy(s, rand(1:s.bits)),                                nothing)
Benchmarks.@benchmarkable(pauliz!,           nothing, pauliz(s, rand(1:s.bits)),                                nothing)
Benchmarks.@benchmarkable(measure!,          nothing, measure(s),                                               nothing)
Benchmarks.@benchmarkable(partialmeasure!,   nothing, partialmeasure(s, rand(1:s.bits), rand(1:s.bits),
                                                                        rand(1:s.bits), rand(1:s.bits),
                                                                        rand(1:s.bits), rand(1:s.bits)),        nothing)
Benchmarks.@benchmarkable(unitarycnot!,      nothing, unitary(s, cn, 4, 3),                                     nothing)
Benchmarks.@benchmarkable(qft!,              nothing, qft(s),                                                   nothing)


benchmarks = [hadamard!, not!, cnot!, ccnot!, swap!, cswap!, sqrtswap!, phaseshift!, paulix!,
              pauliy!, pauliz!, measure!, partialmeasure!, unitarycnot!, qft!]

for b in benchmarks
    println("Benchmark: ", b)
    result = Benchmarks.execute(b)
    println(result)
end
