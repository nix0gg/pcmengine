from physics_constants import G, g
import math

def universal_gravitational_force(m1, m2, r):
    return G * (m1 * m2) / r**2

def gravitational_acceleration(m, r):
    return G * m / r**2

def escape_velocity(m, r):
    return math.sqrt(2 * G * m / r)

def orbital_velocity(m, r):
    return math.sqrt(G * m / r)

def time_period(m, r):
    return 2 * math.pi * math.sqrt(r**3 / (G * m))

def weight(m):
    return m * g

def gravitational_potential_energy(m1, m2, r):
    return -G * (m1 * m2) / r

def gravitational_constant(F,r,m1,m2):
    return F * r**2 / (m1 * m2)

def relationship_between_g_and_G(M,R):
    return 4/3 * math.pi * G * M / R**2

def acceleration_due_to_gravity_on_planet(G,M,R):
    return G * M / R**2

def due_to_altitude(g0,h,R):
    return g0 * (R / (R + h))**2

def due_to_depth(g0,d,R):
    return g0 * (1 - d / R)
__all__ = [name for name in globals() if not name.startswith("_")]
