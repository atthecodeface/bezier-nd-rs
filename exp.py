#!/usr/bin/env python3
import math


def i32_30(x):
    x = round(x * (1 << 30))
    if x >= (1 << 31):
        raise Exception("Number too large for i32_30 fixed point")
    if -x >= (1 << 31):
        raise Exception("Number too negative for i32_30 fixed point")
    return x


def apply_tables(table, e_value):
    exp = 1
    log = 0
    for l, power in table:
        new_exp = exp + (exp / 2 ** (power))
        if e_value > new_exp:
            exp = new_exp
            log += l

    return log


def apply_inv_tables(table, l_value):
    exp = 1
    log = 0
    for l, power in table:
        new_log = log + l
        if l_value > new_log:
            exp = exp + (exp / 2 ** (power))
            log = new_log

    return exp


def apply_inv_tables_2(table, l_value):
    exp = 1
    log = 0
    for l, power in table:
        for x in range(2):
            new_log = log - l
            if new_log <= l_value:
                exp = exp * (2**power - 1) / 2**power
                log = new_log

    return exp


def ln(x):
    return apply_tables(log_table, x)


def exp(x):
    return apply_inv_tables(log_table, x)


def exp_m(x):
    return apply_inv_tables_2(inv_log_table, x)


powers = [i for i in range(31)]

log_table = [(math.log(1 + math.pow(2, -x)), x) for x in powers]
inv_log_table = [
    (math.log((math.pow(2, x + 1) - 1) / math.pow(2, x + 1)), x + 1) for x in powers
]

N = 1000
for i in range(N):
    v = 1 + (i / N)
    log = math.log(v)
    l = ln(v)
    e = exp(l)
    e_m = exp_m(l)
    err_log = math.fabs(log - l)
    err_exp = math.fabs(v - e)
    err_exp_m = math.fabs(e * e_m - 1)
    print(f"{v} {l} {log} {err_log} {err_exp} {err_exp_m} {e} {e_m} {e * e_m}")

    if err_log > 1e-6:
        raise Exception("Log(x) error too large")
    if err_exp > 1e-6:
        raise Exception("Exp(x) error too large")
    if err_exp_m > 1e-6:
        raise Exception("Exp(-x) error too large")
    pass
