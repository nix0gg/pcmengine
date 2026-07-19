import math

def is_reflexive(relation, domain):
    return all((a,a) in relation for a in domain)

def is_symmetric(relation):
    return all((b,a) in relation for (a,b) in relation)

def is_transitive(relation):
    return all((a,c) in relation for (a,b) in relation for (c,d) in relation if b==c)

def is_equivalence(relation,domain):
    return is_reflexive(relation,domain) and is_symmetric(relation) and is_transitive(relation)

def domain(relation):
    return set(a for (a,b) in relation)

def codomain(relation):
    return set(b for (a,b) in relation)

def range_of(relation):
    return codomain(relation)

def inverse_relation(relation):
    return set((b,a) for (a,b) in relation)

def composition(r1,r2):
    return set((a,c) for (a,b) in r1 for (b2,c) in r2 if b == b2)

def is_function(relation,domain):
    return all(sum(1 for (a,b) in relation if a == x) ==1 for x in domain)

def is_injective(relation):
    values = [b for (a,b) in relation]
    return len(values) == len(set(values))

def is_surjective(relation, codomain_set):
    return set((b for (a,b) in relation) == codomain_set)

def is_bijective(relation,domain,codomain_set):
    return is_function(relation,domain) and is_injective(relation) and is_surjective(relation, codomain_set)

def identity_relation(domain):
    return set((a,a) for a in domain)

def number_of_relations(n):
    return 2 ** (n*n)

def number_of_functions(m, n):
    return n**m

def number_of_injections(m,n):
    if m>n:
        return 0
    esult = 1
    for i in range(m):
        result *=(n-i)
    return result

def number_of_bijections(n):
   return math.factorial(n)

__all__ = [name for name in globals() if not name.startswith("_")]