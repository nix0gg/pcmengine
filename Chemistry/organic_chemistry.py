def retention_factor(distance_travelled_compund,distance_travelled_solvent):
    return distance_travelled_compund / distance_travelled_solvent
global compound_mass_taken
#Qualitative Analysis
def qa_c(co2mass,compound_mass_taken):
    return (12/44) * (co2mass/compound_mass_taken) * 100
    
def qa_h(h20mass,):
    return (2/18)*(h20mass/compound_mass_taken)*100

def qa_n(voln2,):
    return (28/22400)* (voln2/compound_mass_taken)*100

def qa_n_kjeldahl(molarity_acid,vol_acid_used,basicity_acid,):
    return (1.4*molarity_acid*vol_acid_used*basicity_acid)/compound_mass_taken

def qa_x(at_mass_x,mass_agx_formed):
    return (at_mass_x/(108+at_mass_x))*(mass_agx_formed/compound_mass_taken)*100

def qa_s(mass_baso4_formed):
    return 