from physics_constants import G, g
import math

def work(f,d,theta):
    return f*d*math.cos(theta)
global theta
def kinetic_energy(m,v):
    return 0.5*m*(v**2)

def potential_energy(m,h):
    return m*g*h

def total_mechanical_energy(ke,pe):
    return ke+pe

def power(w,t):
    return w/t

def power_2(f,v,theta):
    return f*v*math.cos(theta)

def work_energy_theorem(k2,k1):
    delta_k = k2-k1
    return (delta_k)

def work_2(f,d):
    return f*d
def efficiency(wo,wi):
    return (wo/wi)*100

def gravitational_potential_energy(m,r1,r2,m2):
    return -G*(m*m2)*((1/r2)-(1/r1))

def spring_potential_energy(k,x):
    return 0.5*k*(x**2)

def elastic_potential_energy(k,x):
    return 0.5*k*(x**2)

def potential_energy_spring(k,x):
    return 0.5*k*(x**2)

def kinetic_energy_spring(k,x,xmax):
    return 0.5*k*(xmax**2-x**2)

__all__ = [name for name in globals() if not name.startswith("_")]
