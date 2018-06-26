module BitOps

export swapbit, iszerobit, setbit, getbit,
       clearbit, flipbit, removebit, todecimal

# Guarantee a 64-bit unsigned integer
const ONE = 0x00000001

function todecimal(bindigits::Vector{Int})
    result = 0
    offset = length(bindigits) - 1
    for d in bindigits
        result |= (d << offset)
        offset -= 1
    end
    return result
end

function swapbit(num::Int, from::Int, to::Int)
    temp = xor((num >>> from), (num >>> to)) & ONE
    return xor(num, ((temp << from) | (temp << to)))
end

function iszerobit(num::Int, bit::Int)
    return (num & (ONE << bit)) == 0
end

function getbit(num::Int, bit::Int)
    return (num & (ONE << bit)) >>> bit
end

function setbit(num::Int, bit::Int)
    return num | (ONE << bit)
end

function flipbit(num::Int, bit::Int)
    return xor(num, (ONE << bit))
end

function clearbit(num::Int, bit::Int)
    return num & ~(ONE << bit)
end

function removebit(num::Int, bit::Int)
    upper = (num & (0xffffffff << (bit + 1))) >>> 1
    lower = num & ((ONE << bit) - 1)
    return upper | lower
end

end
