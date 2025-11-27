from Physics.physics_constants import N_A
# I don't know at this point lmfao
def concentration(amt_solute, vol_solution):
    return (amt_solute / vol_solution)
def at_mass(mass_one_atom):
    return mass_one_atom / amu
def gram_atom(mass_element_g,gam):
    return mass_element_g / gam

def rmm(avg_mass):
    return avg_mass / ((1/12)*c12mass)

def no_of_moles(substance_mass, molar_mass):
    return substance_mass / molar_mass

def petit_law(specific_heat):
    return 6.4 / specific_heat

def molecular_mass(avg_mass_mol):
    return avg_mass_mol /((1/12)*c12mass)

def equivalent_mass(atomic_mass, valency):
    return atomic_mass / valency

def equivalent_mass_salt(formula_mass,total_charge):
    return formula_mass / total_charge

def equivalent_mass_acid(molecular_mass, basicity):
    return molecular_mass / basicity

# Mole concept
def no_of_moles1(mass,molecular_mass):
    return mass / molecular_mass

def no_of_moles2(no_of_particles):
    return no_of_particles / N_A

def no_of_moles3(gram_atom_element):
    return  gram_atom_element / N_A

def no_of_moles4(gmmass):
    return gmmass / N_A

def no_of_moles5(volume_gas_litres):
    return (volume_gas_litres / 22.4) * N_A

def percentage_composition(mass,molecular_mass):
    return (mass / molecular_mass) * 100

def simplewholenoratio(molecular_mass,empirical_formula_mass):
    return molecular_mass / empirical_formula_mass

#Stoichiometry
def mass_percentage_solute(mass_solute,mass_solution):
    return (mass_solute / mass_solution) * 100

def mole_fraction(moles_A,moles_B):
    return moles_A / (moles_A + moles_B)

def molarity(no_of_moles,volume_solution_litres):
    return no_of_moles / volume_solution_litres

def molarity2(wb,mb,v):
    return wb / (mb * v)
 
def molality(no_of_moles,weight_solvent_kg):
    return no_of_moles / weight_solvent_kg

def molality2(wb,mb,wa):
    return wb / (mb * wa)

def normality(gram_equivalent,volume):
    return gram_equivalent / volume







































#Constants
amu = 1.66e-24
c12mass = 6.023e23