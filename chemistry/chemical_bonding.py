def formal_charge(valence_electrons, non_bonding_electrons, bonding_electrons):
    f_charge = valence_electrons - non_bonding_electrons -((1/2)*bonding_electrons)
    return f_charge

