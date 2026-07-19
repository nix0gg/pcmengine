import math

def degree_of_unsaturation(C, H, N=0, X=0, O=0):
    return (2 * C + 2 + N - H - X) / 2


def h_count_saturated(C, N=0):
    return 2 * C + 2 + N

def co2_from_combustion(moles_hydrocarbon, C_count):
    return moles_hydrocarbon * C_count  

def h2o_from_combustion(moles_hydrocarbon, H_count):
    return moles_hydrocarbon * (H_count / 2)  

def o2_required_combustion(C, H, O=0):
    return C + H / 4 - O / 2

def heat_of_combustion(moles, delta_H_combustion):
    return moles * delta_H_combustion

def halogen_substitution_products(C_count, H_count, halogen="Cl"):
    return f"C{C_count}H{H_count - 1}{halogen}"

def bromine_water_test(moles_alkene):
    return moles_alkene  

def h2_addition(moles_unsaturated, degree_unsat):
    return moles_unsaturated * degree_unsat

def hx_addition_markovnikov(alkene, HX):
    return f"H adds to carbon with more H; {HX[1:]} adds to carbon with fewer H"

def ozonolysis_products(double_bond_position, C_count):
    return f"Cleaves at C{double_bond_position}-C{double_bond_position + 1} double bond"

def electrophilic_substitution(benzene_moles, reagent):
    reactions = {
        "nitration": "C6H5NO2 (nitrobenzene)",
        "sulphonation": "C6H5SO3H (benzenesulphonic acid)",
        "halogenation": "C6H5X (halobenzene)",
        "alkylation": "C6H5R (alkylbenzene)",
        "acylation": "C6H5COR (aryl ketone)"
    }
    return reactions.get(reagent.lower(), "Unknown reaction")

def baeyers_test(moles_alkene):
    return moles_alkene  


def boiling_point_trend(carbon_chain_length):
    return 20 * carbon_chain_length 

def branching_effect_on_bp(n_isomer_bp, branch_count):
    return n_isomer_bp - (branch_count * 7)


def structural_isomers_butane():
    return ["n-butane: CH3CH2CH2CH3", "isobutane: (CH3)3CH"]

def structural_isomers_pentane():
    return [
        "n-pentane: CH3(CH2)3CH3",
        "isopentane: (CH3)2CHCH2CH3",
        "neopentane: C(CH3)4"
    ]

def has_geometric_isomerism(groups_on_c1, groups_on_c2):
    return len(set(groups_on_c1)) == 2 and len(set(groups_on_c2)) == 2

def acidic_order():
    return "Terminal alkynes > alkenes > alkanes (sp > sp2 > sp3)"

__all__ = [name for name in globals() if not name.startswith("_")]
