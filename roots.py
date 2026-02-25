#!/usr/bin/env python3
import math


def recip_sqrt_i32(x, one):
    # An = 1 / sqrt(x)
    # An^2.x = 1
    # B = A^2.x
    # C = A.x
    # A <= A + (1<<i)
    # C <= C + (x<<i)
    # B <= B + 2.((Ax)<<i) + (x<<2i) = B + C<<(i+1) + (x<<2i)
    #
    # If new_B <=1 then (A,B,C) <= new
    A = 0
    B = 0
    C = 0
    for i in range(30, 0, -1):
        # print(A, B, C, x, i, one)
        new_B = B + ((C << (i + 1)) >> 60) + (x >> (60 - 2 * i))
        new_C = C + (x << i)
        new_A = A + (1 << i)
        # print(new_A, new_B, new_C)
        if new_B <= one:
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
    result = recip_sqrt_i32((v << 30) // 8, 1 << 30)
    print(
        v / 8, math.sqrt(8 / v), result / (1 << 30), (result / (1 << 30)) ** 2 * v / 8
    )
    pass
