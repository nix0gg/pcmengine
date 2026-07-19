import math
# Alum: KAl(SO4)2·12H2O
alum_molar_mass = 474

def formal_charge(valence_e, lone_pair_e, bond_e):
    return valence_e - lone_pair_e - (bond_e / 2)

def electronegativity_similarity(en1, en2):
    return abs(en1 - en2)

def borax_bead_test_colour(metal_oxide):
    colours = {
        "CuO": "blue (oxidising), grey (reducing)",
        "CoO": "blue (both flames)",
        "NiO": "brown",
        "MnO2": "violet",
        "Cr2O3": "green",
        "Fe2O3": "yellow/brown"
    }
    return colours.get(metal_oxide, "Unknown")


def thermite_iron_produced(moles_al):
    return moles_al * 56  

def thermite_alumina_produced(moles_al):

    return (moles_al / 2) * 102  

def alum_al_percent():
    return (27 / alum_molar_mass) * 100

def alum_water_percent():
    return (216 / alum_molar_mass) * 100

allotropes_carbon = ["diamond", "graphite", "fullerene", "graphene"]

def hybridisation(allotrope):
    hyb = {
        "diamond": "sp3",
        "graphite": "sp2",
        "fullerene": "sp2",
        "graphene": "sp2"
    }
    return hyb.get(allotrope.lower(), "Unknown")

def co2_produced_from_carbon(moles_c):

    return moles_c * 44  

def co_produced_from_carbon(moles_c):

    return moles_c * 28  

def sio2_from_silicon(moles_si):

    return moles_si * 60  

def si_o_ratio(silicate_type):
    ratios = {
        "orthosilicate": "1:4",
        "pyrosilicate": "2:7",
        "cyclic": "1:3",
        "single chain": "1:3",
        "double chain": "4:11",
        "sheet": "2:5",
        "framework": "1:2"
    }
    return ratios.get(silicate_type.lower(), "Unknown")

def catenation_order():
    return ["C >> Si > Ge > Sn"]

def max_oxidation_state(group_number):
    return group_number - 10  

def inert_pair_effect_state(group_number):
    return group_number - 12  


def bond_angle_hydride(element):
    angles = {
        "CH4": 109.5,
        "SiH4": 109.5,
        "GeH4": 109.5,
        "SnH4": 109.5,
        "NH3": 107,
        "PH3": 93.6,
        "H2O": 104.5,
    }
    return angles.get(element, "Unknown")

__all__ = [name for name in globals() if not name.startswith("_")]






