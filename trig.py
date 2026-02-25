#!/usr/bin/env python3
import math
from functools import reduce


def apply_tables(table, value, d1, d0):
    for base, delta in table:
        scale = 1
        old_value = value
        if value > 0:  # base / 2:
            scale = -1
            pass
        value += base * scale
        new_d0 = d0 + scale * (d1 / 2 ** (-delta))
        new_d1 = d1 - scale * (d0 / 2 ** (-delta))
        d0 = new_d0
        d1 = new_d1
        # print(".  ", value, d0, d1, old_value, base, math.degrees(base))
        pass
    return (value, (d1, d0))


def sincos(x):
    return apply_tables(sincos_table, x, 0, cos_scale)


def apply_tables_inv(table, d1, d0):
    value = 0
    for atan, power in table:
        for x in range(2):
            scale = 1
            if d1 > 0:
                scale = -1
                pass
            value += atan * scale
            new_d0 = d0 + scale * (d1 / 2 ** (-power))
            new_d1 = d1 - scale * (d0 / 2 ** (-power))
            d0 = new_d0
            d1 = new_d1
            pass

    print(f"Returning {value} {d1} {d0}")
    return (value, (d1, d0))


def atan2(y, x):
    return apply_tables_inv(sincos_table, x, y)


def asin(s):
    c = math.sqrt(1 - s * s)
    return atan2(s, c)


def acos(c):
    s = math.sqrt(1 - c * c)
    return atan2(s, c)


def i32_30(x):
    x = round(x * (1 << 30))
    if x >= (1 << 31):
        raise Exception("Number too large for i32_30 fixed point")
    if -x >= (1 << 31):
        raise Exception("Number too negative for i32_30 fixed point")
    return x


angles = [-i for i in range(31)]
dbl_angles = [-i for j in range(31) for i in [j, j, j, j, j]]
print(dbl_angles)
sincos_table = [(math.atan(math.pow(2, x)), x) for x in angles]
asincos_table = [(math.atan(math.pow(2, x)), x) for x in dbl_angles]

print(f"PI/2 as value 0x{i32_30(math.pi / 2):x}")
print("pub const NEG_POW2_I32_30 :&[u8]= &[", end="")
for t, x in sincos_table:
    t = i32_30(t)
    if t == 0:
        break
    print(f"{-x},")
    pass
print("];")

cos_scale = 1
print("pub const ATAN_ANGLES_I32_30 : &[i32]= &[", end="")
for angle, x in sincos_table:
    angle_i32 = i32_30(angle)
    if angle_i32 == 0:
        break
    cos_scale *= math.cos(angle)
    print(f"0x{angle_i32:x},")
    pass
print("];")
print(f"pub const COS_SCALE_I32_30:i32 = 0x{i32_30(cos_scale):x};")

print(cos_scale, sincos_table)


for angle in range(100):
    angle = math.pi * (angle / 100 / 2)
    sc = sincos(angle)
    remainder = sc[0]
    sin = sc[1][0]
    cos = sc[1][1]
    d = sin * sin + cos * cos
    err_sin = math.fabs(sin - math.sin(angle))
    err_cos = math.fabs(cos - math.cos(angle))
    print(f"{angle * 180 / math.pi} sin,cos: {sin},{cos} {err_sin},{err_cos}")

    test_angle = -(math.pi / 2 + atan2(sin, cos)[0])
    print(f"{angle * 180 / math.pi} {test_angle} {test_angle * 180 / math.pi}")
    pass
