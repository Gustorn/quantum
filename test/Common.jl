import Base.isapprox
import FactCheck.roughly

using QSpice.State

roughly(q1::QuantumState) = (q2::QuantumState) -> isapprox(q1, q2)
roughly(q1::Vector{Complex{Float64}}) = (q2::QuantumState) -> isapprox(q1, q2.vector)

function isapprox(q1::Vector{Complex{Float64}}, q2::Vector{Complex{Float64}}, atol=0.0001)
    for (x, y) in zip(q1, q2)
        if !isapprox(x, y, atol = atol)
            return false
        end
    end
    return true
end

