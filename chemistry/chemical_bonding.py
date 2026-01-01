def formal_charge(valence_electrons, non_bonding_electrons, bonding_electrons):
    return valence_electrons - (non_bonding_electrons + bonding_electrons // 2)

