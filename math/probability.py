import math

def classical_probability(favourable, total):
    return favourable / total

def complement(p):
    return 1 - p

def probability_not_a(p_a):
    return 1 - p_a

def odds_in_favour(p):
    return p / (1 - p)

def odds_against(p):
    return (1 - p) / p

def probability_from_odds_favour(m, n):

    return m / (m + n)

def probability_from_odds_against(m, n):
    return n / (m + n)

def p_a_or_b(p_a, p_b, p_a_and_b):
    return p_a + p_b - p_a_and_b

def p_a_or_b_mutually_exclusive(p_a, p_b):
    return p_a + p_b

def p_a_and_b_independent(p_a, p_b):
    return p_a * p_b


def conditional_probability(p_a_and_b, p_b):
    return p_a_and_b / p_b

def p_a_and_b_from_conditional(p_a_given_b, p_b):
    return p_a_given_b * p_b


def bayes(p_b_given_a, p_a, p_b):

    return (p_b_given_a * p_a) / p_b

def total_probability(p_bi_list, p_a_given_bi_list):
    return sum(p * q for p, q in zip(p_bi_list, p_a_given_bi_list))


def ncr(n, r):
    return math.factorial(n) // (math.factorial(r) * math.factorial(n - r))

def prob_exactly_r_from_n(n, r, p):

    q = 1 - p
    return ncr(n, r) * (p ** r) * (q ** (n - r))

def prob_at_least_one(n, p):
    return 1 - (1 - p) ** n

def expected_value_binomial(n, p):
    return n * p

def variance_binomial(n, p):
    return n * p * (1 - p)

def std_dev_binomial(n, p):
    return math.sqrt(variance_binomial(n, p))

def sample_space_coins(n):
    return 2 ** n

def sample_space_dice(n):
    return 6 ** n

def sample_space_cards():
    return 52

__all__ = [name for name in globals() if not name.startswith("_")]
