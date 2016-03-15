module Util

export swap_bit, is_zero, set_bit,
       clear_bit, flip_bit,
       bit_range_equal, shift_range_down

# Executes the XOR swap on two individual bits
# `from` and `to` are the 0 based indicies of the bit
function swap_bit(num::Int, from::Int, to::Int)
    # Store the XOR'd value of the two bits we want to swap
    # in the LSB of temp
    temp = ((num >>> from) $ (num >>> to)) & 0x01
    # First shift the bits we want to swap in the correct position,
    # then XOR them with the original number to complete the swap
    return num $ ((temp << from) | (temp << to))
end

# Checks if a specific bit in a number is zero.
# `bit` is the 0-based index of the bit we want to check
function is_zero(num::Int, bit::Int)
    return (num & (0x01 << bit)) == 0
end

# Sets the bit indicated by `bit` to 1 and returns the new value
function set_bit(num::Int, bit::Int)
    return num | (0x01 << bit)
end

function flip_bit(num::Int, bit::Int)
    return num $ (0x01 << bit)
end

# Sets the bit indicated by `bit` to 0 and returns the new value
function clear_bit(num::Int, bit::Int)
    return num & ~(0x01 << bit)
end

function bit_range_equal(num::Int, start::Int, pattern::Int)
    return ((num & (pattern << start)) >>> start) == pattern
end

# Constructs the number that we get from removing the specified
# bit from the original number, and shifting the upper half down
# 1, while keeping the lower half intact
function shift_range_down(num::Int, bit::Int)
    upper = ((num & (0xffffffff << (bit + 1))) >>> 1)
    lower = num & ((0x01 << bit) - 1)
    return upper | lower
end


end
