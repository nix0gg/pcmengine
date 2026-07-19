import math

def binomial_nCr(n, r):
    return math.factorial(n) // (math.factorial(r) * math.factorial(n - r))

def binomial_term(n, r, a, b):
    return binomial_nCr(n, r) * (a ** (n - r)) * (b ** r)

def binomial_expansion(n, a, b):
    return [binomial_term(n, r, a, b) for r in range(n + 1)]

def general_term(n, r, a, b):
    return binomial_term(n, r, a, b)

def middle_term(n, a, b):
    if n % 2 == 0:
        r = n // 2
        return (f"T_{r+1}", binomial_term(n, r, a, b))
    else:
        r1, r2 = (n - 1) // 2, (n + 1) // 2
        return (f"T_{r1+1} and T_{r2+1}",
                binomial_term(n, r1, a, b),
                binomial_term(n, r2, a, b))

def coefficient_of_xn(n, a, b, power, x_power_in_a=1, x_power_in_b=1):
    for r in range(n + 1):
        x_exp = (n - r) * x_power_in_a + r * x_power_in_b
        if x_exp == power:
            return binomial_term(n, r, a, b)
    return 0

def term_independent_of_x(n, a, b, x_power_in_a, x_power_in_b):
    denom = x_power_in_a - x_power_in_b
    if denom == 0:
        return None
    r = (n * x_power_in_a) / denom
    if r == int(r) and 0 <= r <= n:
        return binomial_term(n, int(r), a, b)
    return None

def binomial_sum(n):
    return 2 ** n

def binomial_alternating_sum(n):
    return 0

def binomial_sum_odd_even(n):
    return 2 ** (n - 1)

def pascal_row(n):
    return [binomial_nCr(n, r) for r in range(n + 1)]

def greatest_coefficient_term(n):
    r = n // 2
    return r, binomial_nCr(n, r)

def numerically_greatest_term(n, a, b, x_val):
    terms = [abs(binomial_term(n, r, a, b)) for r in range(n + 1)]
    r = terms.index(max(terms))
    return r + 1, terms[r]

__all__ = [name for name in globals() if not name.startswith("_")]
