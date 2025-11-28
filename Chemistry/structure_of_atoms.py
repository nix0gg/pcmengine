import math
from Physics.physics_constants import c, e

e_me_ratio = 1.758820*10^11

def frequency(lbda):
    return c / lbda

def velocity(lbda, f):
    return f * lbda

def wave_number(f, c):
    return f / c

#Line spectra
def lyman_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 1**2) - (1 / n2**2))

def balmer_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 2**2) - (1 / n2**2))

def paschen_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 3**2) - (1 / n2**2))

def brackett_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 4**2) - (1 / n2**2))

def pfund_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 5**2) - (1 / n2**2))

def humphreys_series(n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / 6**2) - (1 / n2**2))

def energy_incident_light(h,lbda):
    return (h * c) / lbda

#Bohr's model derivations

def bohr_energy(n,z):
    return (-1312*(z**2))/(n**2)

def bohr_radius(n,z):
    return (0.529*(n**2))/z

def bohr_speed(n,z):
    return (2.18*10**6 * z)/n

def line_spectra_hydrogen(n1,n2):
    R = 1.097373e7  # m^-1
    return R * ((1 / n1**2) - (1 / n2**2))
def debroglie_equation1(h,p):
    return h / p

def debroglie_equation2(h,m,v):
    return h / (m * v)

def debroglie_equation3(h,m,q,v):
    return h / math.sqrt(2*m*q*v)

__all__ = [name for name in globals() if not name.startswith("_")]
