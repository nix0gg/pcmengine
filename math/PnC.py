import math
from collections import Counter 

def factorial(n):
    return math.factorial(n)

def nPr(n,r):
    if r > n:
        return 0
    else:
        return math.factorial(n) // (math.factorial(n-r))
    
def nCr(n,r):
    if r>n:
        return 0
    else:
        return math.factorial(n) // (math.factorial(r)*math.factorial(n-r))

def nCr_relation(n,r):
    return(n,n-r)

def nPr_from_nCr(n,r):
    return nCr(n,r)*math.factorial(r)

def permutation_repetition(n,r):
    return n**r

def combination_repition(n,r):
    return nCr(n+r-1,r)

def permutation_circular(n):
    return math.factorial(n-1)

def circular_permutation_necklace(n):
    return math.factorial(n-1)//2

def permutation_objectsareidentical(n,*groups):
    denom = 1
    for g in groups:
        denom*=math.factorial(g)
    return math.factorial(n)//denom

def selection_minimumis1(n):
    return 2**n-1

def totalsubsets(n):
    return 2**n

def divisorcount(n):
    count = 0
    for i in range(1, int(math.sqrt(n)+1)):
        if n % i == 0:
            if i !=n // i:
                count+=2
            else:
                count = 1

    return count

def word_letter(letters):
    freq = Counter(letters)
    n = len(letters)
    denom = 1
    for f in freq.values():
        denom*= math.factorial(f)
    return math.factorial(n) //denom

def nCr_sum_row(n):
    return 2**n

__all__ = [name for name in globals() if not name.startswith("_")]