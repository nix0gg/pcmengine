altered = {
    '<': '>', 
    '>': '<', 
    '<=': '>=', 
    '>=': '<=',
    }


def solve_linear_inequality(a, b, c, symbol):
    rhs = (c - b) / a
    if a < 0:
        
        symbol = altered.get(symbol, symbol)
    if symbol == '>':
        return f"x > {rhs}"
    elif symbol == '<':
        return f"x < {rhs}"
    elif symbol == '>=':
        return f"x >= {rhs}"
    elif symbol == '<=':
        return f"x <= {rhs}"

def solution_set_intersection(l1, u1, l2, u2):
  
    lower = max(l1, l2)
    upper = min(u1, u2)
    if lower > upper:
        return None 
    return (lower, upper)

def solution_set_union(l1, u1, l2, u2):
    return (min(l1, l2), max(u1, u2))

def is_in_solution(x, lower, upper, inclusive_lower=True, inclusive_upper=True):
    left = (x >= lower) if inclusive_lower else (x > lower)
    right = (x <= upper) if inclusive_upper else (x < upper)
    return left and right

def double_inequality(a, expr_val, b):
    return a < expr_val < b

def solve_absolute_inequality_lt(a, c):
    return (a - c, a + c)

def solve_absolute_inequality_gt(a, c):
    return (a - c, a + c)  

def number_of_integer_solutions(lower, upper):
    import math
    return max(0, math.floor(upper) - math.ceil(lower) + 1)

__all__ = [name for name in globals() if not name.startswith("_")]
