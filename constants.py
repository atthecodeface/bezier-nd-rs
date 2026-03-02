#!/usr/bin/env python3
import math
from cmath import pi


def num_bits(n):
    bit = 1
    while n > 0:
        bit += 1
        n >>= 1
        pass
    return bit


def sqrt(x_sq):
    n = num_bits(x_sq)

    # An = sqrt(x)
    # An^2 = x
    # B = A^2 (int.(2*frac_bits))
    # A <= A + (1<<i)
    # B <= B + 2.(A<<i) + (1<<2i) = B + A<<(i+1) + (1<<2i)
    #
    # If new_B <=1 then (A,B,C) <= new
    A = 0
    B = 0

    for i in range(n, -1, -1):
        # print(A, B, C, x, i, one)
        new_B = B + (A << (i + 1)) + (1 << (2 * i))
        new_A = A + (1 << i)
        if new_B <= x_sq:
            B = new_B
            A = new_A
            pass
        pass
    return A


def recip_sqrt_i32(x, frac_bits):

    # An = 1 / sqrt(x)
    # An^2.x = 1
    # B = A^2.x << (3*fract_bits)
    # C = A.x
    # A <= A + (1<<i)
    # C <= C + (x<<i)
    # B <= B + 2.((Ax)<<i) + (x<<2i) = B + C<<(i+1) + (x<<2i)
    #
    # If new_B <=1 then (A,B,C) <= new
    A = 0
    B = 0
    C = 0
    for i in range(31):
        # i = 30 - i
        # print(A, B, C, x, i, one)
        new_B = B + (C << (32 - i)) + (x << (62 - 2 * i))
        new_C = C + (x << (31 - i))
        new_A = A + (1 << (31 - i))
        # print(new_A, new_B, new_C)
        if new_B <= (1 << (3 * frac_bits)):
            B = new_B
            A = new_A
            C = new_C
            pass
        pass
    return A


# This uses the BBP digit extraction algorithm from Wikipedia (Simon Plouffe, 1995)
def pi_scaled_by(scale):
    sum = 0
    for k in range(scale // 4):
        s = (
            (4 << scale) // (8 * k + 1)
            - (2 << scale) // (8 * k + 4)
            - (1 << scale) // (8 * k + 5)
            - (1 << scale) // (8 * k + 6)
        )
        sum += s >> (4 * k)
        pass
    return sum


def e_scaled_by_scale(scale):
    sum = 0
    factorial = 1
    k = 0
    while factorial < (16 << scale):
        sum += (1 << scale) // factorial
        factorial *= k + 1
        k += 1
        pass
    return sum


def ln_half_scaled_by_scale(scale):
    # Use the Taylor expansion at 1 with x=-1/2 yields -2^-1 - 2^-2/2 - 2^-3/3 - 2^-4/4 ...
    sum = 0
    for k in range(1, scale):
        sum -= ((1 << scale) // k) >> k
        pass
    return sum


def ln_four_fifths_scaled_by_scale(scale):
    # Use the Taylor expansion at 1 with x=-1/5 yields -5^-1 - 5^-2/2 - 5^-3/3 - 2^-4/4 ...
    sum = 0
    fifths = 5
    k = 1
    while fifths < (5 << scale):
        sum -= (1 << scale) // (fifths * k)
        fifths *= 5
        k += 1
        pass
    return sum


ln_two_2048 = -ln_half_scaled_by_scale(2048)
ln_two_256 = ln_two_2048 >> (2048 - 256)

ln_five_quarters_2048 = -ln_four_fifths_scaled_by_scale(2048)
ln_ten_2048 = ln_five_quarters_2048 + 3 * ln_two_2048
ln_ten_256 = ln_ten_2048 >> (2048 - 256)

log2_e_2048 = (1 << 4096) // ln_two_2048
log2_e_256 = log2_e_2048 >> (2048 - 256)

log10_e_2048 = (1 << 4096) // ln_ten_2048
log10_e_256 = log10_e_2048 >> (2048 - 256)

log2_ten_2048 = (log2_e_2048 << 2048) // log10_e_2048
log2_ten_256 = log2_ten_2048 >> (2048 - 256)

log10_two_2048 = (log10_e_2048 << 2048) // log2_e_2048
log10_two_256 = log10_two_2048 >> (2048 - 256)

diff_ln_ten = math.log(10.0) - ln_ten_256 / (1 << 256)
diff_ln_two = math.log(2.0) - ln_two_256 / (1 << 256)

diff_log2_e = math.log(math.e, 2) - log2_e_256 / (1 << 256)
diff_log10_e = math.log(math.e, 10) - log10_e_256 / (1 << 256)

diff_log2_ten = math.log(10, 2) - log2_ten_256 / (1 << 256)
diff_log10_two = math.log(2, 10) - log10_two_256 / (1 << 256)

e_2048 = e_scaled_by_scale(2048)
e_256 = e_2048 >> (2048 - 256)
diff_e = math.e - e_256 / (1 << 256)

print(f"e - f64.e = {diff_e}")

print(f"ln(2) - f64.ln(2) = {diff_ln_two}")
print(f"ln(10) - f64.ln(10) = {diff_ln_ten}")

print(f"log2(e) - f64.log2(e) = {diff_log2_e}")
print(f"log2(10) - f64.log2(10) = {diff_log2_ten}")

print(f"log10(e) - f64.log10(e) = {diff_log10_e}")
print(f"log10(10) - f64.log2(10) = {diff_log2_ten}")

# print(f"{e_256 * 4:x}")

# PI with 256 fractional bits = 0x3243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C89
pi_2048 = pi_scaled_by(2048)
pi_256 = pi_2048 >> (2048 - 256)
pi_sqrt_1024 = sqrt(pi_2048)
pi_sqrt_128 = (pi_sqrt_1024 >> (1024 - 128)) + 1  # Rounding up is better
two_2048 = 2 << 2048
two_256 = 2 << 256
two_sqrt_1024 = sqrt(two_2048)
two_sqrt_128 = two_sqrt_1024 >> (1024 - 128)
r_pi_2048 = (1 << 4096) // pi_2048
r_pi_256 = r_pi_2048 >> (2048 - 256)
r_pi_sqrt_1024 = sqrt(r_pi_2048)
r_pi_sqrt_128 = r_pi_sqrt_1024 >> (1024 - 128)

diff_pi = (pi_256 - pi_sqrt_128 * pi_sqrt_128) / (1 << 128)
print(f"PI - sqrt(PI)^2 = {diff_pi} * 2^-128")
diff_two = (two_256 - two_sqrt_128 * two_sqrt_128) / (1 << 128)
print(f"Two - sqrt(two)^2 = {diff_two} * 2^-128")
diff_pi_r_pi = (pi_256 * r_pi_256) / (1 << 512)
print(f"PI / r_pi = 1+- {diff_pi_r_pi} * 2^-128. (pi*r_pi = {r_pi_256 * pi_256:x}")

print(f"pub const PI:u128 = 0x{pi_256 >> 130:x}; // PI as u128_126")
print(f"pub const FRAC_PI_3:u128 = 0x{(pi_256 // 3) >> 129:x}; // PI/3 as u128_127")
print(f"pub const FRAC_2_PI:u128 = 0x{r_pi_256 >> 127:x}; // 2/PI as u128_128")
print(f"pub const SQRT_PI:u128 = 0x{pi_sqrt_128 >> 1:x}; // PI.sqrt() as u128_127")
print(
    f"pub const FRAC_1_SQRT_PI:u128 = 0x{r_pi_sqrt_128:x}; // 1/pi.sqrt() as u128_128"
)

print(f"pub const SQRT_2:u128 = 0x{two_sqrt_128 >> 1:x}; // 2.sqrt() as u128_127")

print(
    f"pub const E:u128 = 0x{(e_256 >> 130) + 1:x}; // E as u128_126"
)  # Rounded up is better
print(f"pub const LN_2:u128 = 0x{(ln_two_256 >> 128):x}; // LN_2 as u128_128")
print(
    f"pub const LN10:u128 = 0x{(ln_ten_256 >> 130) + 1:x}; // LN_10 as u128_126"
)  # Rounded up is better

print(
    f"pub const LOG2_E:u128 = 0x{(log2_e_256 >> 129) + 1:x}; // LOG2_E as u128_128"
)  # Rounded up is better
print(
    f"pub const LOG2_10:u128 = 0x{(log2_ten_256 >> 130) + 1:x}; // LOG2_10 as u128_126"
)  # Rounded up is better

print(f"pub const LOG10_E:u128 = 0x{(log10_e_256 >> 127):x}; // LOG10_E as u128_128")
print(f"pub const LOG10_2:u128 = 0x{(log10_two_256 >> 127):x}; // LOG10_2 as u128_126")
