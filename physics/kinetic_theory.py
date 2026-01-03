import math
from physics_constants import kB, R_inf

#Pressure:
def pressure1(n,T,V):
    P = n*((R_inf*T)/V)
    return P

def pressure2(T,V):
    P = (R_inf*T) /V
    return P

#Kinetic Interpretation of Pressure
def kip1(row,vbar):
    P = 1/3*(row*(vbar**2))
    return P

def kip2(P,row):
    vrms = math.sqrt((3*P)/row)
    return vrms

def kip3(m,N,V):
    P = 1/3*((m*N)/V)
    return P

def kip4(m,N,vbar):
    PV = 1/3*(m*N*(vbar**2))
    return PV

def kip5(m,N,vbar):
    PV = 2/3*1/2*(m*N*(vbar**2))
    return PV

def kip6(ke):
    PV = 2/3*ke
    return PV
#Kinetic interpretation of Temperature
def kit(T):
    E = 3/2*R_inf*T
    return E

def abs_0_temp(E):
    PV = 2/3*E
    return PV

def monoatomic_ke(T):
    E = 3/2*kB*T
    return E

def diatomic_ke(T):
    E = 5/2*kB*T
    return E

def internal_energy_monoatomic1(N,T):
    U = 3/2*N*kB*T
    return U

def internal_energy_diatomic2(N,T):
    U = 5/2*N*kB*T
    return U

def internal_energy_monoatomic3(P,V):
    U = 3/2*(P*V)
    return U

def internal_energy_diatomic3(P,V):
    U = 5/2*(P*V)
    return U

#specific heat capacity
def cv():
    Cv = 5/2*R_inf
    return Cv

def cp():
    Cp = 7/2*R_inf
    return Cp

