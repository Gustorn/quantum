module State

using Iterators.imap

import Base: convert, copy, getindex, length,
             print, show, showcompact,
             start, done, next

export QuantumState, QUBIT_0, QUBIT_1,
       from_bloch, from_states, from_vector,
       to_bases, print_bases

# A n-qubit quantum state is represented by the Kronecker-product of their
# qubits, stored in `vector`, and the number of bits that make up that state,
# stored in `bits`
type QuantumState
    vector::Vector{Complex{Float64}}
    bits::Int
end

# The two constants representing the |0> and |1> states (respectively)
const QUBIT_0 = QuantumState([1,0], 1)
const QUBIT_1 = QuantumState([0,1], 1)

# Constructs a new quantum state from the given state vector. It automatically
# calculates the number of bits
function from_vector{T<:Number}(state::Vector{T})
    num_bits = floor(Int, log2(length(state)))
    return QuantumState(state, num_bits)
end

# Constructs a new quantum state by taking the Kronecker product of any number
# of other quantum states
function from_states(state0::QuantumState, states::QuantumState...)
    state = mapfoldl(x -> x.vector, kron, state0.vector, states)
    return from_vector(state)
end

# Returns the 1-qubit quantum state that represents the (r, theta, phi) location
# on the Bloch-sphere
function from_bloch(r::Float64, theta::Float64, phi::Float64)
    state = [cos(angle.theta / 2), exp(im * angle.phi) * sin(angle.theta / 2)]
    return QuantumState(state, 1)
end

# Converts the quantum state into the tuples of its bases and their magnitudes
function to_bases(state::QuantumState)
    n_bits = floor(Int, log2(length(state.vector) - 1))
    to_bitstring = x -> (bits(x[1] - 1)[end - n_bits:end], x[2])
    return imap(to_bitstring, enumerate(state.vector)) |> collect
end

convert(state::QuantumState) = state.vector
copy(state::QuantumState) = QuantumState(copy(state.vector), copy(state.bits))

length(state::QuantumState) = length(state.vector)
getindex(state::QuantumState, i::Int) = state.vector[i]

start(state::QuantumState) = start(state.vector)
done(state::QuantumState, current) = done(state.vector, current)
next(state::QuantumState, current) = next(state.vector, current)

showcompact(state::QuantumState) = showcompact(state.vector)

function print(state::QuantumState)
    println("Bits: ", state.bits)
    print("State: [", state.vector[1])

    for i in 2:length(state.vector)
        print(",", state.vector[i])
    end
    print("]")
end

function print_bases(state::QuantumState)
    bases = to_bases(state)
    longest = maximum(imap(x -> length(string(x[2])), bases))
    for basis in bases
        println(rpad(string(basis[2]), longest),
                "    |", basis[1], ">    ",
                "Probability: ", real(conj(basis[2]) * basis[2]))
    end
end

end
