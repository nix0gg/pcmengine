import math
def delta_internal_energy(w,q):
    return q + w

def q_heat_absorbed(delta_U,w):
    return delta_U - w

def work_done_isothermal(n,R,p1,p2,T1):
    return n * R * T1 * math.log(p1 / p2)

def  work_done_ireversible_isothermal(P_ext,vf,vi):
    return P_ext * (vf - vi)
 
def free_expansion(P_ext,dV):
    return P_ext * dV

def work_done_isothermal_reversible(n,R,V1,V2,T1):
    return n * R * T1 * math.log(V2 / V1)

def enthalpy_change(delta_U,p,delta_V):
    return delta_U + p * delta_V

def relation_dH_dU(dU,dNg,R,T):
    return dU + R * T * dNg

def Dng(n2,n1):
    return n2-n1

def work_done_chemical_rections(P,dV):
    return (-P)*dV

