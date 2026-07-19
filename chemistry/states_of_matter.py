import math

# Boyle's Law: P1V1 = P2V2
def boyles_law_p2(p1, v1, v2):
    return (p1 * v1) / v2

def boyles_law_v2(p1, v1, p2):
    return (p1 * v1) / p2

# Charles' Law: V1/T1 = V2/T2
def charles_law_v2(v1, t1, t2):
    return (v1 * t2) / t1

def charles_law_t2(v1, t1, v2):
    return (v2 * t1) / v1

# Gay-Lussac's Law: P1/T1 = P2/T2
def gay_lussac_p2(p1, t1, t2):
    return (p1 * t2) / t1

def gay_lussac_t2(p1, t1, p2):
    return (p2 * t1) / p1

# Combined Gas Law
def combined_gas_law_v2(p1, v1, t1, p2, t2):
    return (p1 * v1 * t2) / (p2 * t1)

# Ideal Gas Law: PV = nRT
R = 8.314  

def ideal_gas_pressure(n, T, V):
    return (n * R * T) / V

def ideal_gas_volume(n, T, P):
    return (n * R * T) / P

def ideal_gas_moles(P, V, T):
    return (P * V) / (R * T)

def ideal_gas_temperature(P, V, n):
    return (P * V) / (n * R)

# Dalton's Law of Partial Pressures
def partial_pressure(mole_fraction, total_pressure):
    return mole_fraction * total_pressure

def total_pressure(*partial_pressures):
    return sum(partial_pressures)

def mole_fraction_from_pressure(partial_pressure, total_pressure):
    return partial_pressure / total_pressure

# Graham's Law of Diffusion
def rate_diffusion(M1, M2, r2):
    # r1/r2 = sqrt(M2/M1)
    return r2 * math.sqrt(M2 / M1)

def molar_mass_from_diffusion(r1, r2, M1):
    return M1 * (r2 / r1) ** 2

def rate_effusion_ratio(M1, M2):
    return math.sqrt(M2 / M1)

# Kinetic Molecular Theory
def rms_speed(M, T):
    return math.sqrt(3 * R * T / M)

def average_speed(M, T):
    return math.sqrt(8 * R * T / (math.pi * M))

def most_probable_speed(M, T):
    return math.sqrt(2 * R * T / M)

def average_kinetic_energy(T):
    kB = 1.380649e-23
    return (3/2) * kB * T

def kinetic_energy_molar(T):
    return (3/2) * R * T

# van der Waals equation
def vdw_pressure(n, V, T, a, b):
    return ((n * R * T) / (V - n * b)) - (a * n**2 / V**2)

def vdw_volume_correction(n, b):
    return n * b

def vdw_pressure_correction(n, V, a):
    return a * n**2 / V**2

# Compressibility factor
def compressibility_factor(P, V, n, T):
    return (P * V) / (n * R * T)

# Vapour pressure and boiling
def relative_lowering_vp(mole_fraction_solute):
    return mole_fraction_solute

def boiling_point_elevation(kb, m):
    return kb * m

def freezing_point_depression(kf, m):
    return kf * m

# Surface tension and viscosity (conceptual solvers)
def surface_tension(force, length):
    return force / length

def viscosity_coefficient(force, area, velocity_gradient):
    return (force / area) / velocity_gradient

__all__ = [name for name in globals() if not name.startswith("_")]
