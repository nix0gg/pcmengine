import math

def limit_sin_x_over_x():
    return 1

def limit_tan_x_over_x():
    return 1

def limit_1_plus_1_over_n_n():
    return math.e

def limit_exp_x_minus_1_over_x():
    return 1

def limit_log_1_plus_x_over_x():
    return 1

def limit_a_pow_x_minus_1_over_x(a):
    return math.log(a)

def derivative_by_first_principle(f, x, h=1e-7):
    return (f(x + h) - f(x)) / h


def d_xn(n, x):
    return n * x ** (n - 1)

def d_constant(c):
    return 0

def d_sin(x):
    return math.cos(x)

def d_cos(x):
    return -math.sin(x)

def d_tan(x):
    return 1 / math.cos(x) ** 2

def d_cosec(x):
    return -1 / (math.sin(x) * math.tan(x))

def d_sec(x):
    return math.tan(x) / math.cos(x)

def d_cot(x):
    return -1 / math.sin(x) ** 2

def d_ex(x):
    return math.exp(x)

def d_ax(a, x):
    return (a ** x) * math.log(a)

def d_ln(x):
    return 1 / x

def d_log_a(a, x):
    return 1 / (x * math.log(a))


def product_rule(f, g, df, dg):

    return df * g + f * dg

def quotient_rule(f, g, df, dg):

    return (df * g - f * dg) / g ** 2

def chain_rule(df_du, du_dx):
    return df_du * du_dx

def sum_rule(df, dg):
    return df + dg

def scalar_multiple_rule(c, df):
    return c * df


def slope_tangent(f, x):
    return derivative_by_first_principle(f, x)

def slope_normal(f, x):
    m = slope_tangent(f, x)
    return -1 / m

def tangent_line_y(f, df_at_x, x0, x):
    return f(x0) + df_at_x * (x - x0)

def normal_line_y(f, df_at_x, x0, x):
    m_normal = -1 / df_at_x
    return f(x0) + m_normal * (x - x0)


def nth_derivative_xn(n, r):

    if r > n:
        return 0
    result = 1
    for i in range(r):
        result *= (n - i)
    return result

def nth_derivative_sin(n, x):
    return math.sin(x + n * math.pi / 2)

def nth_derivative_cos(n, x):
    return math.cos(x + n * math.pi / 2)

__all__ = [name for name in globals() if not name.startswith("_")]
