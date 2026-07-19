import math
protium_mass = 1
deuterium_mass = 2
tritium_mass = 3


def density_water_at_temp(T_celsius):
    return 1000 * (1 - ((T_celsius - 3.98) ** 2) / (503.57 * (T_celsius + 283)))

def ph_water(T_kelvin):

    kw = 1e-14 * (10 ** ((13.83 - 4787.3 / T_kelvin)))
    return -0.5 * math.log10(kw)


def temporary_hardness_removed(ca_hco3_conc):
    return ca_hco3_conc

def ppm_hardness(mass_calcium_salt_g, volume_water_L, molar_mass_salt, molar_mass_caco3=100):
    moles = mass_calcium_salt_g / molar_mass_salt
    mass_caco3_equivalent = moles * molar_mass_caco3
    return (mass_caco3_equivalent / (volume_water_L * 1000)) * 1e6

def degree_of_hardness_ppm(mass_caco3_equiv_g, volume_water_L):
    return (mass_caco3_equiv_g / (volume_water_L * 1000)) * 1e6

def lime_required(hardness_ppm, volume_L):
    moles_caco3_equiv = (hardness_ppm * volume_L * 1000 * 1e-6) / 100
    mass_lime = moles_caco3_equiv * 74
    return mass_lime


def h2o2_volume_strength(molarity):
    return molarity * 11.2

def h2o2_molarity_from_volume_strength(volume_strength):
    return volume_strength / 11.2

def h2o2_normality(molarity):

    return molarity * 2

def h2o2_molarity_from_normality(normality):
    return normality / 2

def h2o2_mass_percent(molarity, density_g_per_mL):
    mass_per_litre = molarity * 34  
    total_mass_per_litre = density_g_per_mL * 1000
    return (mass_per_litre / total_mass_per_litre) * 100

def volume_h2_produced(moles_water):
    return moles_water * 22.4

def volume_o2_produced(moles_water):
    return (moles_water / 2) * 22.4

def bond_energy_order(bond_energies_dict):
    return sorted(bond_energies_dict.items(), key=lambda x: x[1], reverse=True)

__all__ = [name for name in globals() if not name.startswith("_")]
