def center_of_mass(masses, positions):
    numerator = sum(masses[i] * positions[i] for i in range(len(masses)))
    denominator = sum(masses)
    return numerator / denominator

def linear_momentum(mass, velocity):
    return mass * velocity

def total_momentum(masses, velocities):
    return sum(masses[i] * velocities[i] for i in range(len(masses)))

def moment_of_inertia(masses, distances):
    return sum(masses[i] * distances[i]**2 for i in range(len(masses)))

def torque(i, alpha):
    return i * alpha

def angular_momentum(i, omega):
    return i * omega

def rotational_kinetic_energy(i, omega):
    return 0.5 * i * omega**2

def parallel_axis_theorem(I_cm, m, d):
    return I_cm + m * d**2


#Value of moment of inertia based on shape

def i_ring(m,r,axis):
    if axis == "perpendicular":
        return m * r**2
    elif axis == "diameter":
        return m*(r**2)/2
    else:
        raise ValueError("Invalid axis. Use 'perpendicular' or 'diameter'.")
    
def i_rod(m,l):
    return m*(l**2)/12

def i_circular_disc(m,r,axis):
    if axis == "perpendicular":
        return m*(r**2)/2
    elif axis == "diameter":
        return m*(r**2)/4
    else:
        raise ValueError("Invalid axis. Use 'perpendicular' or 'diameter'.")
    
def i_cylinder(m,r,kind):
    if kind == "solid":
        return m*(r**2)/2
    elif kind == "hollow":
        return m*(r**2)
    else:
        raise ValueError("Invalid kind. Use 'solid' or 'hollow'.")
    
def i_sphere(m,r):
    return 2*(m*(r**2))/5

__all__ = [name for name in globals() if not name.startswith("_")]
