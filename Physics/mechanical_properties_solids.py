import math

#Stress
def longitudinal_stress(F, A): 
    return F / A

def hooke_law_stress(E, strain):
    return E * strain

def shearing_stress1(F, A):
    return F / A

def shearing_stress2(eta, theta):
    return eta * theta

#Strain
def shearing_strain(dx,L):
    return dx / L

def longitudinal_strain(dL, L):
    return dL / L

def volumetric_strain(dV, V):
    return dV / V

#Young's Modulus
def youngs_modulus(stress, strain):
    return stress / strain

def elongation(longitudinal_strain, L, stress, strain):
    return longitudinal_strain * L / youngs_modulus(stress, strain)

def shear_modulus(f, a, dx, l):
    return (f / a) / (dx / l)

def shear_modulus2(f,l,a,dx):
    return (f * l) / (a * dx)

def shear_modulus3(F,A,theta):
    return F / A * math.tan(theta)

def bulk_modulus_elasticity1(F, A, dV, V):
    return ((-F) / A) / (dV / V)

def bulk_modulus_elasticity2(P, dV, V):
    return (-P * V) / dV

def compressibility1(B):
    return 1 / B

def compressibility2(dV, V, P):
    return -(dV / (P * V))

def poisson_ratio(lateral_strain, longitudinal_strain):
    return -lateral_strain / longitudinal_strain

def poissonratio2(dr,r,dl,l):
    return -(dr / r) / (dl / l)

def elastic_potential_energy1(F,dL):
    return (F/2)* dL 

def elastic_potential_energy2(F,A,dL,L):
    return 0.5 * (F/A) * (dL/L) * (A*L)
 
def y1(B,sigma):
    return 3*B *(1- 2*sigma)

def y2(eta,sigma):
    return 2*eta*(1 + sigma)

def y3(sigma,b,eta):
    return (3*b-2*eta)/ (2*eta+6*b)
def y9(B,eta):
    return (1/B) + (3/eta)

def bulk_modulus(stress,strain):
    return stress /strain

def bulk_modulus2(f,a,dv,v):
    return (-f*v) / (a*dv)

def depression(W,L,Y,b,d):
    return (W*L) / (4*Y*b*d**3)


__all__ = [name for name in globals() if not name.startswith("_")]
