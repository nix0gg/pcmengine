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

def enthalpy_change1(delta_U,p,delta_V):
    return delta_U + p * delta_V

def relation_dH_dU(dU,dNg,R,T):
    return dU + R * T * dNg

def Dng(n2,n1):
    return n2-n1

def work_done_chemical_rections(P,dV):
    return (-P)*dV

def heat_capacity(q,dT):
    return q /dT

def specific_heat_capacity(C,m):
    return C / m

def molar_heat_capacity(q,n,dT):
    return q / (n * dT)

def heat_capacity_at_constant_volume(delta_U,dT):
    return delta_U / dT

def heat_capacity_at_constant_pressure(delta_H,dT):
    return delta_H / dT

def heat_combustion(W,m,t2,t1,c,w1,M):
    return (((W+m)*(t2-t1)*c)/w1)*M

def entropy(qrev, T):
    return qrev / T

def gibbs_energy(H,T,S):
    return H - T*S

def gibbs_helmholtz_energy(dH,T,dS):
    return dH - T*dS

def gibbs_function(dG0, R, T, K):
     return -2.303 * R * T * math.log(K)

def temperature_equilibrium(dH,dS):
    return dH / dS

def equilibrium(dH,T,dS):
    return dH - T*dS

def isochoric_change(cv,dT):
    return cv*dT

def enthalpy_change2(dU,p,dV):
    return dU + p*dV

def spontaneous_process(dS_sys,dS_surr):
    return dS_sys + dS_surr

