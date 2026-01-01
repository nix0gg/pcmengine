import math

def viscosity(force, area, velocity_gradient):
    return (force / area) / velocity_gradient

def bernoulli_theorem(P,row,g,h,v):
    return P + row * g * h + 0.5 * row * v**2

def torricelli_law(g,h,p,p0,row):
    return math.sqrt(2 * g * h * (p - p0) / row)

def excess_pressure1(S,R):
    return 2*S / R

def excess_pressure2(S,R):
    return 4*S /R

def gauge_pressure(row,g,h):
    return row * g * h

def hydraulic_lift(a2,a1,f1):
    return (a2 / a1) * f1

def terminal_velocity(r,row,row0,g,eta):
    return (2 * r**2 * g * (row0 - row)) / (9 * eta)

def u_tube(S,row,g,r1,r2):
    return (2*S / (row * g))  * (r1 - r2)

def stoke_law(eta,r,v):
    return 6 * math.pi * eta * r * v

def capillary_rise(S, row, g, r, theta):
    return 2*S*math.cos(theta) / (row * g * r)

def absolute_pressure(P0,h,row,g):
    return P0 + row * g * h

__all__ = [name for name in globals() if not name.startswith("_")]
