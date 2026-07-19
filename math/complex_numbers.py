import math

def complex_add(a, b, c, d):
    return (a + c, b + d)

def complex_subtract(a, b, c, d):
    return (a - c, b - d)

def complex_multiply(a, b, c, d):
    return (a * c - b * d, a * d + b * c)
    

def complex_divide(a, b, c, d):
    denom = c ** 2 + d ** 2
    return ((a * c + b * d) / denom, (b * c - a * d) / denom)

def complex_modulus(a, b):
    return math.sqrt(a ** 2 + b ** 2)

def argument(a, b):
    return math.degrees(math.atan2(b, a))

def conjugate(a, b):
    return (a, -b)

def polar_form(a, b):
    r = complex_modulus(a, b)
    theta = argument(a, b)
    return (r, theta)

def rectangular_form(r, theta_deg):
    theta = math.radians(theta_deg)
    return (r * math.cos(theta), r * math.sin(theta))

def power(a, b, n):
    r = complex_modulus(a, b)
    theta = argument(a, b)
    rn = r ** n
    theta_n = math.radians(n * theta)
    return (rn * math.cos(theta_n), rn * math.sin(theta_n))

def nth_root(a, b, n, k=0):
    r = complex_modulus(a, b)
    theta = math.radians(argument(a, b))
    r_root = r ** (1 / n)
    theta_k = (theta + 2 * math.pi * k) / n
    return (r_root * math.cos(theta_k), r_root * math.sin(theta_k))

def complex_modulus_squared(a, b):
    return a ** 2 + b ** 2

def multiplicative_inverse(a, b):
    mod_sq = complex_modulus_squared(a, b)
    return (a / mod_sq, -b / mod_sq)

def sum_of_conjugates(a, b):
    return 2 * a

def product_with_conjugate(a, b):
    return complex_modulus_squared(a, b)

def discriminant(a, b, c):
    return b ** 2 - 4 * a * c

def quadratic_roots_complex(a, b, c):
    d = discriminant(a, b, c)
    if d >= 0:
        r1 = (-b + math.sqrt(d)) / (2 * a)
        r2 = (-b - math.sqrt(d)) / (2 * a)
        return ((r1, 0), (r2, 0))
    else:
        real_part = -b / (2 * a)
        imag_part = math.sqrt(-d) / (2 * a)
        return ((real_part, imag_part), (real_part, -imag_part))

__all__ = [name for name in globals() if not name.startswith("_")]