#!/usr/bin/env python3
import math


def sqrt_i32(x, frac_bits):
    # An = sqrt(x)
    # An^2 = x
    # B = A^2 (int.(2*frac_bits))
    # A <= A + (1<<i)
    # B <= B + 2.(A<<i) + (1<<2i) = B + A<<(i+1) + (1<<2i)
    #
    # If new_B <=1 then (A,B,C) <= new
    A = 0
    B = 0

    for i in range(31):
        # print(A, B, C, x, i, one)
        new_B = B + (A << (32 - i)) + (1 << (62 - 2 * i))
        new_A = A + (1 << (31 - i))
        # print(new_A, new_B, new_C)
        if new_B <= x << frac_bits:
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


def i32_30(x):
    x = round(x * (1 << 30))
    if x >= (1 << 31):
        raise Exception("Number too large for i32_30 fixed point")
    if -x >= (1 << 31):
        raise Exception("Number too negative for i32_30 fixed point")
    return x


N = 100
for x in range(N):
    v = 17 * (x + 1)
    result = recip_sqrt_i32((v << 30) // 8, 30)
    result2 = sqrt_i32(v << 23, 22)
    print(
        2 * v,
        v / 8,
        math.sqrt(8 / v),
        result / (1 << 30),
        (result / (1 << 30)) ** 2 * v / 8,
        (result2 * result2) / (1 << 44),
    )
    err = (result / (1 << 30)) ** 2 * v / 8 - 1
    if abs(err) > 1e-6:
        raise Exception("Error in reciprocal square root")
    pass
