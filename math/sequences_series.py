import math

def ap_xterm(a, d, n):
    return a + (n - 1) * d

def ap_sum(a, d, n):
    return n * (2 * a + (n - 1) * d) / 2

def ap_sum_from_last(n, a, l):
    return n * (a + l) / 2

def ap_common_difference(a, b):
    return b - a

def ap_xterm_from_last(l, d, n):
    return l - (n - 1) * d

def ap_number_of_terms(a, d, l):
    return int((l - a) / d) + 1

def arithmetic_mean(a, b):
    return (a + b) / 2

def arithmetic_means_between(a, b, n):
    d = (b - a) / (n + 1)
    return [a + i * d for i in range(1, n + 1)]


def gp_xterm_term(a, r, n):
    return a * r ** (n - 1)

def gp_sum(a, r, n):
    if r == 1:
        return a * n
    return a * (r ** n - 1) / (r - 1)

def gp_sum_infinite(a, r):
    if abs(r) >= 1:
        return None  
    return a / (1 - r)

def gp_common_ratio(a, b):
    return b / a

def geometric_mean(a, b):
    return math.sqrt(a * b)

def geometric_means_between(a, b, n):
    r = (b / a) ** (1 / (n + 1))
    return [a * r ** i for i in range(1, n + 1)]

def gp_product_of_terms(a, r, n):

    return (a ** n) * (r ** (n * (n - 1) // 2))


def hp_xterm(a, d, n):

    return 1 / (a + (n - 1) * d)

def harmonicmean(a, b):
    return 2 * a * b / (a + b)


def am_gm_hm(a, b):
    am = arithmetic_mean(a, b)
    gm = geometric_mean(a, b)
    hm = harmonic_mean(a, b)
    return am, gm, hm  

def sum_of_n(n):
    return n * (n + 1) // 2

def sum_of_n_squared(n):
    return n * (n + 1) * (2 * n + 1) // 6

def sum_of_n_cubed(n):
    return (n * (n + 1) // 2) ** 2

def sum_of_squares_diff(n):
    return sum_of_n_squared(n)

def method_of_differences(sequence, n):
    return sum(sequence[:n])

__all__ = [name for name in globals() if not name.startswith("_")]
