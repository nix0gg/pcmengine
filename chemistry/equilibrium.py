import math

# Equilibrium constant expressions
def kc(product_concs, product_coeffs, reactant_concs, reactant_coeffs):
    numerator = 1
    for c, n in zip(product_concs, product_coeffs):
        numerator *= c ** n
    denominator = 1
    for c, n in zip(reactant_concs, reactant_coeffs):
        denominator *= c ** n
    return numerator / denominator

def kp(kc, dng, R, T):
    return kc * (R * T) ** dng

def kc_from_kp(kp, dng, R, T):
    return kp / (R * T) ** dng

def dng(moles_products, moles_reactants):
    return moles_products - moles_reactants

# Degree of dissociation
def degree_of_dissociation(moles_dissociated, initial_moles):
    return moles_dissociated / initial_moles

def degree_of_dissociation_kp(kp, P):
    # For A(g) ⇌ B(g) + C(g) type
    return math.sqrt(kp / (kp + P))

# Relation between Kc and degree of dissociation (for A ⇌ nB)
def kc_from_alpha(alpha, C, n):
    # alpha = degree of dissociation, C = initial conc, n = moles of product
    return (alpha ** n * C ** (n - 1)) / (1 - alpha)

# pH calculations
def ph(H_conc):
    return -math.log10(H_conc)

def poh(OH_conc):
    return -math.log10(OH_conc)

def h_conc_from_ph(pH):
    return 10 ** (-pH)

def oh_conc_from_poh(pOH):
    return 10 ** (-pOH)

def ph_from_poh(pOH):
    return 14 - pOH

def poh_from_ph(pH):
    return 14 - pH

# Ionic product of water
def kw(H_conc, OH_conc):
    return H_conc * OH_conc

# Acid dissociation constant
def ka(H_conc, A_conc, HA_conc):
    return (H_conc * A_conc) / HA_conc

def ka_from_alpha(alpha, C):
    return (alpha ** 2 * C) / (1 - alpha)

def alpha_from_ka(ka, C):
    return math.sqrt(ka / C)

def ph_weak_acid(ka, C):
    return 0.5 * (- math.log10(ka) - math.log10(C))

# Base dissociation constant
def kb(OH_conc, BH_conc, B_conc):
    return (OH_conc * BH_conc) / B_conc

def kb_from_alpha(alpha, C):
    return (alpha ** 2 * C) / (1 - alpha)

def poh_weak_base(kb, C):
    return 0.5 * (-math.log10(kb) - math.log10(C))

# Relation between Ka and Kb
def kw_from_ka_kb(ka, kb):
    return ka * kb

def ka_from_kw_kb(kw, kb):
    return kw / kb

def kb_from_kw_ka(kw, ka):
    return kw / ka

# Buffer solutions (Henderson-Hasselbalch)
def ph_buffer_acid(pka, salt_conc, acid_conc):
    return pka + math.log10(salt_conc / acid_conc)

def ph_buffer_base(pkb, salt_conc, base_conc):
    return 14 - pkb - math.log10(salt_conc / base_conc)

# Solubility product
def ksp_ab(a_conc, b_conc):
    return a_conc * b_conc

def ksp_ab2(a_conc, b_conc):
    return a_conc * b_conc ** 2

def solubility_from_ksp_1_1(ksp):
    return math.sqrt(ksp)

def solubility_from_ksp_1_2(ksp):
    # AB2 type: Ksp = 4s^3
    return (ksp / 4) ** (1/3)

# Reaction quotient
def reaction_quotient(product_concs, product_coeffs, reactant_concs, reactant_coeffs):
    return kc(product_concs, product_coeffs, reactant_concs, reactant_coeffs)

# van't Hoff factor
def vant_hoff_factor(observed_colligative, theoretical_colligative):
    return observed_colligative / theoretical_colligative

__all__ = [name for name in globals() if not name.startswith("_")]
