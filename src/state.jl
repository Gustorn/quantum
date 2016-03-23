module State

using Iterators

import Base: convert, copy, getindex, setindex!, length,
             show, start, done, next, endof, eltype, isapprox

export QuantumState, QUBIT0, QUBIT1, EMPTY_STATE, BELL_STATE,
       frombloch, fromstates, fromvector, randomstate

# A n-qubit quantum state is represented by the Kronecker-product of their
# qubits, stored in `vector`, and the number of bits that make up that state,
# stored in `bits`
type QuantumState
    vector::Vector{Complex{Float64}}
    bits::Int
end

# The two constants representing the |0> and |1> states (respectively)
const QUBIT0 = QuantumState([1,0], 1)
const QUBIT1 = QuantumState([0,1], 1)

# TODO(gustorn): maybe get rid of the empty state and go for Nullable{QuantumState}
# where applicable
const EMPTY_STATE = QuantumState([], 0)
const BELL_STATE  = QuantumState(1 / sqrt(2) .* [1,0,0,1], 2)

# Constructs a new quantum state from the given state vector. It automatically
# calculates the number of bits
function fromvector{T<:Number}(state::Vector{T})
    numbits = floor(Int, log2(length(state)))
    return QuantumState(state, numbits)
end

# Constructs a new quantum state by taking the Kronecker product of any number
# of other quantum states
function fromstates(state0::QuantumState, states::QuantumState...)
    state::Vector{Complex{Float64}} = mapfoldl(x -> x.vector, kron, state0.vector, states)
    return fromvector(state)
end

# Returns the 1-qubit quantum state that represents the (r, theta, phi) location
# on the Bloch-sphere
function frombloch(theta::Float64, phi::Float64)
    state = [cos(theta / 2), exp(im * phi) * sin(theta / 2)]
    return QuantumState(state, 1)
end

function isapprox(q1::QuantumState, q2::QuantumState, atol=0.0001)
    for (x, y) in zip(q1, q2)
        if !isapprox(x, y, atol = atol)
            return false
        end
    end
    return true
end

function randomstate(nbits::Int)
    vec = [rand() + rand() * im for _ in 1:2^nbits]
    n = sqrt(sumabs2(vec))
    vec = vec .* (1.0 / n)
    return QuantumState(vec, nbits)
end

convert(state::QuantumState) = state.vector
copy(state::QuantumState) = QuantumState(copy(state.vector), copy(state.bits))

length(state::QuantumState) = length(state.vector)
getindex(state::QuantumState, i::Int) = state.vector[i]
setindex!(state::QuantumState, c::Complex{Float64}, i::Int) = state.vector[i] = c

start(state::QuantumState) = start(state.vector)
done(state::QuantumState, current) = done(state.vector, current)
next(state::QuantumState, current) = next(state.vector, current)
endof(state::QuantumState) = endof(state.vector)
eltype(::Type{QuantumState}) = Complex{Float64}

function show(io::IO, state::QuantumState)
    for i = 1:length(state)
        basis = i - 1

        real = Base.real(state[i])
        resign = real < 0.0 ? '-' : ' '
        real = abs(real)

        imag = Base.imag(state[i])
        imsign = imag < 0.0 ? '-' : '+'
        imag = abs(imag)

        @printf("%c%04.4f %c %04.4fi    |%s>    P: %04.4f\n",
                resign, real, imsign, imag, bin(basis, state.bits), abs2(state[i]))
    end
end

end
