from physics_constants import N_A, R_inf
#First Law of Thermodynamics
def firstlaw1(dU,dW):
    dQ = dU + dW
    return dQ

def firstlaw2(dQ,dW):
    dU = dQ - dW
    return dU

def firstlaw3(P,dV):
    dW = P*dV
    return dW

def firstlaw4(dU,P,dV):
    dQ = dU + P * dV
    return dQ

#Specific heat capacity
def specific_heat_capacity1(m,s,dT):
    q = m*s*dT
    return q

def speicif_heat_capacity2(m,s,tf,ti):
    q = m*s*(tf-ti)
    return q

def specific_heat_capacity3(q,m,dT):
    s = q / (m*dT)
    return s

def molar_specific_heat_solid1(H,mew):
    C = H / mew
    return C

def molar_specific_heat_solid2(mew,dQ,dT):
    C = (1/mew)*(dQ/dT)
    return C

def average_energy_oscillator_1d(kb,T):
    kbt = (0.5*kb*T)*2
    return kbt

#Total Energy for a mole of a solid
def total_energy_mole_solid1(kb,T):
    U = 3 * kb * T * N_A
    return U

def total_energy_mole_solid2(T):
    U = 3 * R_inf * T
    return U

def total_energy_mole_solid3(dQ,dt):
    C = dQ / dt
    return C

def total_energy_mole_solid4(dU,dt):
    C = dU / dt
    return C

def total_energy_mole_solid():
    C = 3* R_inf
    return C


#Meyer's Relation

def meyereqn1(dQ,dT):
    Cv = dQ / dT
    return Cv

def meyereqn1_5(dU,dT):
    Cv = dU / dT 
    return Cv

def meyereqn2(dQ,dT):
    Cp = dQ/dT
    return Cp

def meyereqn2_1(dU,dT,P,dV):
    Cp = (dU / dT) + P*(dV/dT)
    return Cp

def meyereqn2_2(Cv,P,dV,dT):
    Cp = Cv + P*(dV/dT)
    return Cp

def meyer_eqn2_3(P,dV,dT):
    R = P*(dV / dT)
    return R

#Internal Energy

def internalenergy_dU(n,Cv,dT):
    dU = n*Cv*dT
    return dU

def internalenergy_dQ(n,Cp,dT):
    dQ = n*Cp*dT
    return dQ


















