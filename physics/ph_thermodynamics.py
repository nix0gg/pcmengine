from physics_constants import N_A, R_inf
import math
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

def dU_cyclic_process(Uf,Ui):
    dU = Uf - Ui 
    return dU

#Equation of states (Same formula, just slightly modified):
#P1*V1 = P2*V2

#P1 = P2V2/V1
def eqn_of_state1(p2,v2,v1):
    p1 = (p2*v2) / v1
    return p1

def eqn_of_state2(p1,v1,v2):
    p2 = (p1*v1)/v2
    return p2

def eqn_of_state3(p2,v2,p1):
    v1 = (p2*v2) / p1
    return v1

def eqn_of_state4(p1,v1,p2):
    v2 = (p1*v1) / p2
    return v2

#Specific heat
def specific_heat(Q,m,dT):
    C = Q /(m*dT)
    return C

#Work done (isothermal process)
def work_isothermal1(mew,T,Pi,Pf):
    W1 = mew*R_inf*math.log10(Pi/Pf)
    return W1

def work_isothermal2(mew,T,Vf,Vi):
    W2 = mew*R_inf*T*math.log10(Vf/Vi)
    return W2

def eqn_of_state_adiabatic1(Cp,Cv):
    gamma = Cp / Cv
    return gamma

def eqn_of_state_adiabatic2(gamma,Cp):
    Cv = Cp / gamma
    return Cv

def eqn_of_state_adiabatic3(Cv,gamma):
    Cp = gamma*Cv
    return Cp

def work_done_adiabatic1(Pi,Vi,Pf,Vf,gamma):
    W3 = (Pi*Vi-Pf*Vf)/(gamma-1)
    return W3

def work_done_adiabatic2(mew,Ti,Tf,gamma):
    W4 = (mew*R_inf*(Ti-Tf))/(gamma-1)
    return W4

#P*V**gamma = K
def work_adiabatic1(K,V,gamma):
    P = K / (V**gamma)
    return P


def work_adiabatic2(K,P):
    v_gamma = K / P
    return v_gamma

#Pf*Vf = mew*R*T
def work_adiabatic3(mew,T,Vf):
    Pf = (mew*R_inf*T) / Vf
    return Pf

def work_adiabatic4(Pf,mew,T):
    Vf = (mew*R_inf*T)/Pf
    return Vf

def work_adiabatic5(Pf,Vf,T):
    mew = (Pf*Vf) / (R_inf *T)
    return mew

#Pi*Vi = mew*R*T
def work_adiabatic6(mew,T,Vi):
    Pi = (mew*R_inf*T) / Vi
    return Pi

def work_adiabatic7(Pi,mew,T):
    Vi = (mew*R_inf*T)/Pi
    return Vi

def work_adiabatic8(Pi,Vi,T):
    mew = (Pi*Vi) / (R_inf *T)
    return mew

def work_done_isochoric1(P,Vi,Vf):
    W = P*(Vf-Vi)
    return W

def work_done_isochoric2(mew,Tf,Ti):
    W = mew*R_inf*(Tf-Ti)
    return W

def work_done_isochoric3(mew,dT):
    W = mew*R_inf*dT
    return W

def work_done_isochoric4(mew,dT):
    W = mew*R_inf*dT
    return W

#Carnot's Engine
#Isothermal Expansion
def isothermal_expansion(T1,V2,V1):
    Q1 = R_inf*T1*math.log(V2/V1)
    return Q1
#
def adiabtatic_expansion1(T1,T2,gamma):
    W = ((R_inf*T2)-(R_inf*T1)) / (1-gamma)
    return W

def adiabatic_expansion2(T1,T2,gamma):
    W = (R_inf*(T1-T2))/(gamma-1)
    return W 

#Isothermal compression
def isothermal_compression1(T2,V4,V3):
    W = -(R_inf)*T2*math.log(V4/V3)
    return W

def isothermal_compression2(T2,V3,V4):
    W = R_inf*T2*math.log(V3/V4)
    return W

#Adiabatic compression

def adiabatic_compression(T1,T2,gamma):
    W = R_inf(T1-T2)/(gamma-1)
    return W

#Efficiency of carnot engine
def eta_carnot_engine1(W,Q1):
    eta = W/Q1
    return eta

def eta_carnot_engine2(Q2,Q1):
    eta = 1-(Q2/Q1)
    return eta

def eta_carnot_engine3(T2,T1,V1,V2,V3,V4):
    eta = (R_inf*T2*math.log(V3/V4)) /(R_inf*T1*math.log(V2/V1))
    return eta

def eta_carnot_engine4(T2,T1):
    eta = 1-(T2/T1)
    return eta




