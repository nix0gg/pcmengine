import math

# Ionization energy trend (comparative, not absolute)
def ionization_energy_trend(IE1, IE2):
    return IE2 / IE1

# Electronegativity difference and bond character
def ionic_character_percent(electronegativity_diff):
    return (1 - math.exp(-0.25 * electronegativity_diff ** 2)) * 100

# Hydration enthalpy
def hydration_enthalpy(charge, radius_pm):
    return charge / radius_pm

# Solubility of Group 2 sulphates (comparative)
def lattice_enthalpy_effect(lattice_enthalpy, hydration_enthalpy):
    return hydration_enthalpy - lattice_enthalpy

# Flame test wavelengths (nm) - constants
flame_wavelengths = {
    "Li": 670,   # crimson red
    "Na": 589,   # golden yellow
    "K":  404,   # lilac/violet
    "Rb": 780,   # red
    "Cs": 455,   # blue
    "Ca": 622,   # brick red
    "Sr": 606,   # crimson
    "Ba": 554,   # apple green
}

def flame_colour(element):
    return flame_wavelengths.get(element, "Unknown")

# Electrolysis - Down's process (NaCl)
def nacl_electrolysis_sodium(moles_nacl):
    return moles_nacl * 23 

def nacl_electrolysis_chlorine(moles_nacl):
    return (moles_nacl / 2) * 71 

# Solvay process - Na2CO3 production
def solvay_nahco3(moles_nacl):
    return moles_nacl * 84

def solvay_na2co3(moles_nahco3):

    return (moles_nahco3 / 2) * 106  

def decomposition_temp_trend(cation_radius):
    return cation_radius  

# Diagonal relationship (Li-Mg, Be-Al)
def charge_density(charge, radius_pm):
    return charge / (radius_pm ** 2)

# Lime (CaO) reactions
def lime_water_reaction(moles_co2):
    return moles_co2 * 100  

def quicklime_from_limestone(moles_caco3):
    return moles_caco3 * 56  

def slaked_lime_from_quicklime(moles_cao):
    return moles_cao * 74 

# Plaster of Paris
def pop_water_of_crystallisation(moles_pop):
 return moles_pop * 1.5 

__all__ = [name for name in globals() if not name.startswith("_")]
