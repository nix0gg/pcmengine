import os
import sys

 # If running the file directly (not as a package), make project root importable
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from physics.physics_constants import g
   

# ---------------------------- Basic Quantities ----------------------------
def proc_area():
    from physics.dimensional_formulas import area
    length = float(input("Enter length (L): "))
    breadth = float(input("Enter breadth (L): "))
    print("Area (L^2):", area(length, breadth))

def proc_volume():
    from physics.dimensional_formulas import volume
    length = float(input("Enter length (L): "))
    breadth = float(input("Enter breadth (L): "))
    height = float(input("Enter height (L): "))
    print("Volume (L^3):", volume(length, breadth, height))

def proc_mass_density():
    from physics.dimensional_formulas import mass_density
    mass = float(input("Enter mass (M): "))
    volume = float(input("Enter volume (L^3): "))
    print("Mass Density (M/L^3):", mass_density(mass, volume))

def proc_frequency():
    from physics.dimensional_formulas import frequency
    time_period = float(input("Enter time period (T): "))
    print("Frequency (1/T):", frequency(time_period))

def proc_speed():
    from physics.dimensional_formulas import speed
    distance = float(input("Enter distance (L): "))
    time = float(input("Enter time (T): "))
    print("Speed (L/T):", speed(distance, time))

def proc_acceleration():
    from physics.dimensional_formulas import acceleration
    change_in_velocity = float(input("Enter change in velocity (L/T): "))
    time = float(input("Enter time (T): "))
    print("Acceleration (L/T^2):", acceleration(change_in_velocity, time))

# ---------------------------- Mechanics ----------------------------
def proc_force():
    from physics.dimensional_formulas import force
    mass = float(input("Enter mass (M): "))
    acceleration = float(input("Enter acceleration (L/T^2): "))
    print("Force (MLT^-2):", force(mass, acceleration))

def proc_impulse():
    from physics.dimensional_formulas import impulse
    force = float(input("Enter force (MLT^-2): "))
    time = float(input("Enter time (T): "))
    print("Impulse (MLT^-1):", impulse(force, time))

def proc_work():
    from physics.dimensional_formulas import work
    force = float(input("Enter force (MLT^-2): "))
    distance = float(input("Enter distance (L): "))
    print("Work (ML^2T^-2):", work(force, distance))

def proc_power():
    from physics.dimensional_formulas import power
    work = float(input("Enter work done (ML^2T^-2): "))
    time = float(input("Enter time taken (T): "))
    print("Power (ML^2T^-3):", power(work, time))

def proc_momentum():
    from physics.dimensional_formulas import momentum
    mass = float(input("Enter mass (M): "))
    velocity = float(input("Enter velocity (L/T): "))
    print("Momentum (MLT^-1):", momentum(mass, velocity))

def proc_pressure():
    from physics.dimensional_formulas import pressure
    force = float(input("Enter force (MLT^-2): "))
    area = float(input("Enter area (L^2): "))
    print("Pressure (ML^-1T^-2):", pressure(force, area))

def proc_strain():
    from physics.dimensional_formulas import strain
    change = float(input("Enter change in dimension (L): "))
    original = float(input("Enter original dimension (L): "))
    print("Strain (dimensionless):", strain(change, original))

def proc_modulus_of_elasticity():
    from physics.dimensional_formulas import modulus_of_elasticity
    stress = float(input("Enter stress (ML^-1T^-2): "))
    strain = float(input("Enter strain (dimensionless): "))
    print("Modulus of Elasticity (ML^-1T^-2):", modulus_of_elasticity(stress, strain))

def proc_surface_tension():
    from physics.dimensional_formulas import surface_tension
    force = float(input("Enter force (MLT^-2): "))
    length = float(input("Enter length (L): "))
    print("Surface Tension (MT^-2):", surface_tension(force, length))

def proc_kinetic_energy():
    from physics.dimensional_formulas import kinetic_energy
    mass = float(input("Enter mass (M): "))
    velocity = float(input("Enter velocity (L/T): "))
    print("Kinetic Energy (ML^2T^-2):", kinetic_energy(mass, velocity))

def proc_potential_energy():
    from physics.dimensional_formulas import potential_energy
    mass = float(input("Enter mass (M): "))
    height = float(input("Enter height (L): "))
    print("Potential Energy (ML^2T^-2):", potential_energy(mass, height, g))

def proc_escape_velocity():
    from physics.dimensional_formulas import escape_velocity
    mass = float(input("Enter mass of body (M): "))
    radius = float(input("Enter radius (L): "))
    G = float(input("Enter gravitational constant (L^3M^-1T^-2): "))
    print("Escape Velocity (L/T):", escape_velocity(mass, radius, G))

# ---------------------------- Electricity ----------------------------
def proc_charge_dimensional_formula():
    from physics.dimensional_formulas import charge
    current = float(input("Enter current (A): "))
    time = float(input("Enter time (T): "))
    print("Charge (AT):", charge(current, time))

def proc_current_density():
    from physics.dimensional_formulas import current_density
    current = float(input("Enter current (A): "))
    area = float(input("Enter area (L^2): "))
    print("Current Density (AL^-2):", current_density(current, area))

def proc_voltage():
    from physics.dimensional_formulas import voltage
    work = float(input("Enter work done (ML^2T^-2): "))
    charge = float(input("Enter charge (AT): "))
    print("Voltage (ML^2T^-3A^-1):", voltage(work, charge))

def proc_resistance():
    from physics.dimensional_formulas import resistance
    voltage = float(input("Enter voltage (ML^2T^-3A^-1): "))
    current = float(input("Enter current (A): "))
    print("Resistance (ML^2T^-3A^-2):", resistance(voltage, current))

# ---------------------------- Optics and Radiation ----------------------------
def proc_refractive_index():
    from physics.dimensional_formulas import refractive_index
    speed_in_medium = float(input("Enter speed in medium (L/T): "))
    print("Refractive Index (dimensionless):", refractive_index(speed_in_medium))

def proc_wavenumber():
    from physics.dimensional_formulas import wavenumber
    wavelength = float(input("Enter wavelength (L): "))
    print("Wavenumber (L^-1):", wavenumber(wavelength))

# ---------------------------- Heat and Thermodynamics ----------------------------
def proc_heat_capacity():
    from physics.dimensional_formulas import heat_capacity
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    temp_change = float(input("Enter temperature change (K): "))
    print("Heat Capacity (ML^2T^-2K^-1):", heat_capacity(heat_energy, temp_change))

def proc_specific_heat_capacity():
    from physics.dimensional_formulas import specific_heat_capacity
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    mass = float(input("Enter mass (M): "))
    temp_change = float(input("Enter temperature change (K): "))
    print("Specific Heat Capacity (L^2T^-2K^-1):", 
          specific_heat_capacity(heat_energy, mass, temp_change))

def proc_latent_heat():
    from physics.dimensional_formulas import latent_heat
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    mass = float(input("Enter mass (M): "))
    print("Latent Heat (L^2T^-2):", latent_heat(heat_energy, mass))

# ---------------------------- Motion in One Dimension ----------------------------


def proc_velocity_1d():
    from physics.motion_in_one_dimension import velocity
    u = float(input("Enter initial velocity (L/T): "))
    t = float(input("Enter time (T): "))
    direction = input("Is direction upward? (yes/no): ").lower()
    isDirectionUp = direction == "yes"
    print("Final Velocity (L/T):", velocity(u, t, isDirectionUp))

def proc_displacement_1d():
    from physics.motion_in_one_dimension import displacement
    u = float(input("Enter initial velocity (L/T): "))
    t = float(input("Enter time (T): "))
    direction = input("Is direction upward? (yes/no): ").lower()
    isDirectionUp = direction == "yes"
    print("Displacement (L):", displacement(u, t, isDirectionUp))
    

def proc_max_height():
    from physics.motion_in_one_dimension import max_height
    u = float(input("Enter initial velocity (L/T): "))
    print("Maximum Height (L):", max_height(u))

def proc_time_of_flight_1d():
    from physics.motion_in_one_dimension import time_of_flight
    u = float(input("Enter initial velocity (L/T): "))
    print("Time of Flight (T):", time_of_flight(u))

def proc_ground_strike_velocity():
    from physics.motion_in_one_dimension import ground_strike_velocity
    h = float(input("Enter height (L): "))
    print("Ground Strike Velocity (L/T):", ground_strike_velocity(h))

def proc_average_velocity():
    from physics.motion_in_one_dimension import average_velocity
    u = float(input("Enter initial velocity (L/T): "))
    v = float(input("Enter final velocity (L/T): "))
    print("Average Velocity (L/T):", average_velocity(u, v))

def proc_average_acceleration():
    from physics.motion_in_one_dimension import average_acceleration
    u = float(input("Enter initial velocity (L/T): "))
    v = float(input("Enter final velocity (L/T): "))
    t = float(input("Enter time (T): "))
    print("Average Acceleration (L/T^2):", average_acceleration(u, v, t))

def proc_normal_force():
    from physics.motion_in_one_dimension import n1
    m = float(input("Enter mass (M): "))
    print("Normal Force (MLT^-2):", n1(m))

def proc_normal_force_accelerated():
    from physics.motion_in_one_dimension import n2,n3
    m = float(input("Enter mass (M): "))
    a = float(input("Enter acceleration (L/T^2): "))
    up = input("Is acceleration upward? (yes/no): ").lower()
    if up == "yes":
        print("Normal Force (MLT^-2):", n2(m, a))
    else:
        print("Normal Force (MLT^-2):", n3(m, a))

# ---------------------------- Motion in Plane ----------------------------
def proc_projectile_range():
    from physics.motion_in_plane import range_of_projectile
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Range of Projectile (L):", range_of_projectile(u, theta))

def proc_projectile_max_height():
    from physics.motion_in_plane import max_height
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Maximum Height of Projectile (L):", max_height(u, theta))

def proc_projectile_time_of_flight():
    from physics.motion_in_plane import time_of_flight
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Time of Flight (T):", time_of_flight(u, theta))

def proc_horizontal_displacement():
    from physics.motion_in_plane import horizontal_displacement
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Horizontal Displacement (L):", horizontal_displacement(u, theta, t))

def proc_vertical_displacement():
    from physics.motion_in_plane import vertical_displacement
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Vertical Displacement (L):", vertical_displacement(u, theta, t))

def proc_resultant_velocity():
    from physics.motion_in_plane import resultant_velocity
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Resultant Velocity (L/T):", resultant_velocity(u, theta, t))

def proc_angle_of_velocity():
    from physics.motion_in_plane import angle_of_velocity
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter initial angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Angle of Velocity (degrees):", angle_of_velocity(u, theta, t))

# ---------------------------- Work, Power and Energy ----------------------------


def proc_work_angular():
    from physics.work_power_energy import work
    force = float(input("Enter force (MLT^-2): "))
    distance = float(input("Enter distance (L): "))
    theta = float(input("Enter angle between force and displacement (radians): "))
    print("Work Done (ML^2T^-2):", work(force, distance, theta))

def proc_total_mechanical_energy():
    from physics.work_power_energy import total_mechanical_energy
    ke = float(input("Enter kinetic energy (ML^2T^-2): "))
    pe = float(input("Enter potential energy (ML^2T^-2): "))
    print("Total Mechanical Energy (ML^2T^-2):", total_mechanical_energy(ke, pe))

def proc_power_from_force():
    from physics.work_power_energy import power_2
    force = float(input("Enter force (MLT^-2): "))
    velocity = float(input("Enter velocity (L/T): "))
    theta = float(input("Enter angle between force and velocity (radians): "))
    print("Power (ML^2T^-3):", power_2(force, velocity,theta))

def proc_work_energy_theorem():
    from physics.work_power_energy import work_energy_theorem
    k2 = float(input("Enter final kinetic energy (ML^2T^-2): "))
    k1 = float(input("Enter initial kinetic energy (ML^2T^-2): "))
    print("Change in Kinetic Energy (ML^2T^-2):", work_energy_theorem(k2, k1))

def proc_efficiency():
    from physics.work_power_energy import efficiency
    work_out = float(input("Enter work output (ML^2T^-2): "))
    work_in = float(input("Enter work input (ML^2T^-2): "))
    print("Efficiency (%):", efficiency(work_out, work_in))

def proc_gravitational_potential_energy():
    from physics.work_power_energy import gravitational_potential_energy
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r1 = float(input("Enter initial separation (L): "))
    r2 = float(input("Enter final separation (L): "))
    print("Gravitational Potential Energy (ML^2T^-2):", 
          gravitational_potential_energy(m1, r1, r2, m2))

def proc_spring_potential_energy():
    from physics.work_power_energy import spring_potential_energy
    k = float(input("Enter spring constant (MT^-2): "))
    x = float(input("Enter displacement from equilibrium (L): "))
    print("Spring Potential Energy (ML^2T^-2):", spring_potential_energy(k, x))

def proc_spring_kinetic_energy():
    from physics.work_power_energy import kinetic_energy_spring
    k = float(input("Enter spring constant (MT^-2): "))
    x = float(input("Enter current displacement (L): "))
    xmax = float(input("Enter maximum displacement (L): "))
    print("Spring Kinetic Energy (ML^2T^-2):", kinetic_energy_spring(k, x, xmax))

# ---------------------------- Mechanical Properties of Solids ----------------------------
def proc_longitudinal_stress():
    from physics.mechanical_properties_solids import longitudinal_stress
    force = float(input("Enter force (MLT^-2): "))
    area = float(input("Enter cross-sectional area (L^2): "))
    print("Longitudinal Stress (ML^-1T^-2):", longitudinal_stress(force, area))

def proc_hooke_law_stress():
    from physics.mechanical_properties_solids import hooke_law_stress
    E = float(input("Enter value of 'E':"))
    strain = float(input("Enter value of strain:"))
    return hooke_law_stress(E,strain)

def proc_shearing_stress1():
    from physics.mechanical_properties_solids import shearing_stress1
    F= float(input("Enter force:"))
    A = float(input("Enter area:"))
    return shearing_stress1(F,A) 

def proc_shearing_stress2():
    from physics.mechanical_properties_solids import shearing_stress2
    eta = float(input("Enter eta:"))
    theta = float(input("Enter theta:"))
    return shearing_stress2(eta,theta)

def proc_shearing_strain():
    from physics.mechanical_properties_solids import shearing_strain
    dx = float(input("Enter dx:"))
    L = float(input("Enter l:"))
    return shearing_strain(dx,L)

def proc_strain_longitudinal():
    from physics.mechanical_properties_solids import longitudinal_strain
    dL = float(input("Enter change in length (m): "))
    L = float(input("Enter original length (m): "))
    result = longitudinal_strain(dL, L)
    print(f"Longitudinal Strain = {result:.5g}")

def proc_shear_modulus():
    from physics.mechanical_properties_solids import shear_modulus
    f = float(input("Enter Force (N): "))
    a = float(input("Enter Area (m²): "))
    dx = float(input("Enter displacement (m): "))
    l = float(input("Enter original length (m): "))
    result = shear_modulus(f, a, dx, l)
    print(f"Shear Modulus = {result:.5g} N/m²")

def proc_bulk_modulus_elasticity2():
    from physics.mechanical_properties_solids import bulk_modulus_elasticity2
    P = float(input("Enter Pressure (Pa): "))
    dV = float(input("Enter change in Volume (m³): "))
    V = float(input("Enter original Volume (m³): "))
    result = bulk_modulus_elasticity2(P, dV, V)
    print(f"Bulk Modulus = {result:.5g} N/m²")

def proc_poisson_ratio():
    from physics.mechanical_properties_solids import poisson_ratio
    lateral_strain = float(input("Enter lateral strain: "))
    longitudinal_strain = float(input("Enter longitudinal strain: "))
    result = poisson_ratio(lateral_strain, longitudinal_strain)
    print(f"Poisson's Ratio = {result:.5g}")

def proc_poisson_ratio2():
    from physics.mechanical_properties_solids import poissonratio2
    dr = float(input("Enter value for dr:"))
    r = float(input("Enter value for r:"))
    dl = float(input("Enter value for dl:"))
    l = float(input("Enter value for l:"))
    return poissonratio2(dr,r,dl,l)

def proc_elastic_potential_energy2():
    from physics.mechanical_properties_solids import elastic_potential_energy2
    F = float(input("Enter Force (N): "))
    A = float(input("Enter Area (m²): "))
    dL = float(input("Enter change in length (m): "))
    L = float(input("Enter original length (m): "))
    result = elastic_potential_energy2(F, A, dL, L)
    print(f"Elastic Potential Energy = {result:.5g} J")

def proc_elastic_potential_energy1():
    from physics.mechanical_properties_solids import elastic_potential_energy1
    F = float(input("Enter value for 'F':"))
    dL = float(input("Enter value for dL:"))
    return elastic_potential_energy1(F,dL)

def proc_volumetric_strain():
        from physics.mechanical_properties_solids import volumetric_strain
        dV = float(input("Enter value of dV:"))
        V = float(input("Enter value for V:"))
        return volumetric_strain(dV,V)

def proc_young_modulus():
    from physics.mechanical_properties_solids import youngs_modulus
    stress = float(input("Enter value for stress:"))
    strain = float(input("Enter value for strain:"))
    return youngs_modulus(stress,strain)

def proc_elongation():
    from physics.mechanical_properties_solids import elongation
    longitudinal_strain = float(input("Enter value for longitudinal strain:"))
    L = float(input("Enter value for l:"))
    stress = float(input("Enter value for stress:"))
    strain= float(input("Enter value for strain:"))
    return elongation(longitudinal_strain,L,stress,strain)

def proc_shear_modulus2():
    from physics.mechanical_properties_solids import shear_modulus2
    f = float(input("Enter value for 'f':"))
    l = float(input("Enter value for 'l':"))
    a = float(input("Enter value for 'a':"))
    dx = float(input("Enter value for 'dx':"))
    return shear_modulus2(f,l,a,dx)

def proc_shear_modulus3():
    from physics.mechanical_properties_solids import shear_modulus3
    F = float(input("Enter value of 'F':"))
    A = float(input("Enter value for 'A':"))
    theta = float(input("Enter value for theta:"))
    return shear_modulus3(F,A,theta)

def proc_modulus_elasticity1():
    from physics.mechanical_properties_solids import bulk_modulus_elasticity1
    F = float(input("Enter value for 'F':"))
    A = float(input("Enter value for 'A'"))
    dV = float(input("Enter value for 'dV':"))
    V = float(input("Enter value for 'V':"))
    return bulk_modulus_elasticity1(F,A,dV,V)

def proc_compressibility1():
    from physics.mechanical_properties_solids import compressibility1
    B  = float(input("Enter value for 'B':"))
    return compressibility1(B)

def proc_compressibility2():
    from physics.mechanical_properties_solids import compressibility2
    dV = float(input("Enter the value of 'dV':"))
    V = float(input("Enter the value of 'V':"))
    P = float(input("Enter value for 'P':"))
    return compressibility2(dV,V,P)

def proc_y1():
    from physics.mechanical_properties_solids import y1
    B = float(input("Enter the value of 'B':"))
    sigma = float(input("Enter the value of sigma:"))
    return y1(B,sigma)

def proc_y2():
    from physics.mechanical_properties_solids import y2
    eta = float(input("Enter the value of eta:"))
    sigma = float(input("Enter the value of sigma:"))
    return y2(eta,sigma)

def proc_y3():
    from physics.mechanical_properties_solids import y3
    sigma = float(input("Enter the value for sigma:"))
    b = float(input("Enter the value for 'b':"))
    eta = float(input("Enter the value for eta:"))
    return y3(sigma,b,eta)

def proc_y9():
    from physics.mechanical_properties_solids import y9
    B = float(input("Enter the value of 'B':"))
    eta = float(input("Enter the value of eta:"))
    return y9(B,eta)

def proc_bulk_modulus():
    from physics.mechanical_properties_solids import bulk_modulus
    stress = float(input("Enter the value for stress:"))
    strain = float(input("Enter the value for strain:"))
    return bulk_modulus(stress,strain)

def proc_bulk_modulus2():
    from physics.mechanical_properties_solids import bulk_modulus2
    f = float(input("Enter the value for 'f':"))
    a = float(input("Enter the value for 'a':"))
    dv = float(input("Enter the value for 'dv':"))
    v = float(input("Enter the value for 'v':"))
    return bulk_modulus2(f,a,dv,v)

def proc_depression():
    from physics.mechanical_properties_solids import depression
    W = float(input("Enter the value for 'W':"))
    L = float(input("Enter the value for 'L':"))
    Y = float(input("Enter the value for 'Y':"))
    b = float(input("Enter the value for 'b':"))
    d = float(input("Enter the value for 'd':"))
    return depression(W,L,Y,b,d)

# ---------------------------- Gravitation ----------------------------
def proc_gravitational_force():
    from physics.gravitation import universal_gravitational_force
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r = float(input("Enter distance between centers (L): "))
    print("Gravitational Force (MLT^-2):", universal_gravitational_force(m1, m2, r))

def proc_gravitational_acceleration():
    from physics.gravitation import gravitational_acceleration
    m = float(input("Enter mass of the body (M): "))
    r = float(input("Enter distance from center (L): "))
    print("Gravitational Acceleration (LT^-2):", gravitational_acceleration(m, r))

def proc_orbital_velocity():
    from physics.gravitation import orbital_velocity
    m = float(input("Enter mass of central body (M): "))
    r = float(input("Enter orbital radius (L): "))
    print("Orbital Velocity (LT^-1):", orbital_velocity(m, r))

def proc_orbital_period():
    from physics.gravitation import time_period
    m = float(input("Enter mass of central body (M): "))
    r = float(input("Enter orbital radius (L): "))
    print("Orbital Period (T):", time_period(m, r))

def proc_weight():
    from physics.gravitation import weight
    m = float(input("Enter mass (M): "))
    print("Weight (MLT^-2):", weight(m))

def proc_gravitational_pe():
    from physics.gravitation import gravitational_potential_energy
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r = float(input("Enter distance between centers (L): "))
    print("Gravitational Potential Energy (ML^2T^-2):", gravitational_potential_energy(m1, m2, r))

# ==========================================================
#  MECHANICAL PROPERTIES OF FLUIDS
# ==========================================================

def proc_viscosity():
    from physics.mechanical_properties_fluids import viscosity
    F = float(input("Enter Force (N): "))
    A = float(input("Enter Area (m²): "))
    dv_dy = float(input("Enter velocity gradient (s⁻¹): "))
    result = viscosity(F, A, dv_dy)
    print(f"Viscosity = {result:.5g} Pa·s")

def proc_bernoulli():
    from physics.mechanical_properties_fluids import bernoulli_theorem
    P = float(input("Enter Pressure (Pa): "))
    rho = float(input("Enter Density (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    h = float(input("Enter height (m): "))
    v = float(input("Enter velocity (m/s): "))
    result = bernoulli_theorem(P, rho, g, h, v)
    print(f"Total Bernoulli Constant = {result:.5g} J/m³")

def proc_torricelli():
    from physics.mechanical_properties_fluids import torricelli_law
    g = float(input("Enter g (m/s²): "))
    h = float(input("Enter height (m): "))
    p = float(input("Enter pressure p (Pa): "))
    p0 = float(input("Enter reference pressure p₀ (Pa): "))
    rho = float(input("Enter density (kg/m³): "))
    result = torricelli_law(g, h, p, p0, rho)
    print(f"Velocity (Torricelli's Law) = {result:.5g} m/s")

def proc_excess_pressure():
    from physics.mechanical_properties_fluids import excess_pressure1
    S = float(input("Enter Surface Tension (N/m): "))
    R = float(input("Enter Radius (m): "))
    result = excess_pressure1(S, R)
    print(f"Excess Pressure = {result:.5g} Pa")

def proc_capillary_rise():
    from physics.mechanical_properties_fluids import capillary_rise
    S = float(input("Enter Surface Tension (N/m): "))
    rho = float(input("Enter Density (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    r = float(input("Enter Radius of tube (m): "))
    theta = float(input("Enter contact angle θ (radians): "))
    result = capillary_rise(S, rho, g, r, theta)
    print(f"Capillary Rise = {result:.5g} m")

def proc_terminal_velocity():
    from physics.mechanical_properties_fluids import terminal_velocity
    r = float(input("Enter Radius of sphere (m): "))
    rho = float(input("Enter density of fluid (kg/m³): "))
    rho0 = float(input("Enter density of body (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    eta = float(input("Enter Viscosity (Pa·s): "))
    result = terminal_velocity(r, rho, rho0, g, eta)
    print(f"Terminal Velocity = {result:.5g} m/s")

def proc_hydraulic_lift():
    from physics.mechanical_properties_fluids import hydraulic_lift
    a2 = float(input("Enter Area₂ (m²): "))
    a1 = float(input("Enter Area₁ (m²): "))
    f1 = float(input("Enter Force₁ (N): "))
    result = hydraulic_lift(a2, a1, f1)
    print(f"Force on larger piston = {result:.5g} N")

def proc_absolute_pressure():
    from physics.mechanical_properties_fluids import absolute_pressure
    P0 = float(input("Enter Atmospheric Pressure (Pa): "))
    h = float(input("Enter height (m): "))
    rho = float(input("Enter Density (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    result = absolute_pressure(P0, h, rho, g)
    print(f"Absolute Pressure = {result:.5g} Pa")


# ==========================================================
# THERMAL PROPERTIES OF MATTER
# ==========================================================

def proc_temperature_conversion():
    from physics.thermal_properties_of_matter import temperature_conversion1
    f = float(input("Enter temperature in Fahrenheit (°F): "))
    result = temperature_conversion1(f)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_temperature_from_kelvin():
    from physics.thermal_properties_of_matter import temperature_conversion2
    k = float(input("Enter temperature in Kelvin (K): "))
    result = temperature_conversion2(k)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_temperature_from_rankine():
    from physics.thermal_properties_of_matter import temperature_conversion3
    R = float(input("Enter temperature in Rankine (°R): "))
    result = temperature_conversion3(R)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_thermal_linear_expansion():
    from physics.thermal_properties_of_matter import thermal_linear_expansion1
    l0 = float(input("Enter original length (m): "))
    beta = float(input("Enter coefficient of linear expansion (1/°C): "))
    delta_T = float(input("Enter change in temperature (°C): "))
    result = thermal_linear_expansion1(l0, beta, delta_T)
    print(f"Final length = {result:.5g} m")

def proc_thermal_volume_expansion():
    from physics.thermal_properties_of_matter import thermal_volume_expansion
    V0 = float(input("Enter original volume (m³): "))
    gamma = float(input("Enter coefficient of volume expansion (1/°C): "))
    delta_T = float(input("Enter temperature change (°C): "))
    result = thermal_volume_expansion(V0, gamma, delta_T)
    print(f"Final volume = {result:.5g} m³")

def proc_latent_heat_fusion():
    from physics.thermal_properties_of_matter import latent_heat_fusion
    m = float(input("Enter mass (kg): "))
    lf = float(input("Enter latent heat of fusion (J/kg): "))
    result = latent_heat_fusion(m, lf)
    print(f"Heat absorbed/released = {result:.5g} J")

def proc_specific_heat_capacity_tpm():
    from physics.thermal_properties_of_matter import specific_heat_capacity
    Q = float(input("Enter heat absorbed (J): "))
    m = float(input("Enter mass (kg): "))
    delta_T = float(input("Enter temperature change (°C): "))
    result = specific_heat_capacity(Q, m, delta_T)
    print(f"Specific Heat Capacity = {result:.5g} J/kg·°C")

def proc_rate_of_heat_flow():
    from physics.thermal_properties_of_matter import rate_of_heat_flow1
    k = float(input("Enter thermal conductivity (W/m·K): "))
    A = float(input("Enter Area (m²): "))
    delta_T = float(input("Enter temperature difference (°C): "))
    d = float(input("Enter thickness (m): "))
    result = rate_of_heat_flow1(k, A, delta_T, d)
    print(f"Rate of Heat Flow = {result:.5g} W")

def proc_newton_law_cooling():
    from physics.thermal_properties_of_matter import newton_law_cooling
    k = float(input("Enter cooling constant: "))
    m = float(input("Enter mass (kg): "))
    s = float(input("Enter specific heat capacity (J/kg·°C): "))
    t2 = float(input("Enter initial temperature (°C): "))
    t1 = float(input("Enter final temperature (°C): "))
    result = newton_law_cooling(k, m, s, t2, t1)
    print(f"Rate of Cooling = {result:.5g} J/s")

# ==========================================================
# ELECTRICITY
# ==========================================================

def proc_charge():
    from physics.dimensional_formulas import charge
    current = float(input("Enter current (A): "))
    time = float(input("Enter time (s): "))
    result = charge(current, time)
    print(f"Charge = {result:.5g} C")

def proc_current_density_electricity():
    from physics.dimensional_formulas import current_density
    current = float(input("Enter current (A): "))
    area = float(input("Enter area (m²): "))
    result = current_density(current, area)
    print(f"Current Density = {result:.5g} A/m²")

def proc_voltage_electricity():
    from physics.dimensional_formulas import voltage
    work = float(input("Enter work done (J): "))
    q = float(input("Enter charge (C): "))
    result = voltage(work, q)
    print(f"Voltage = {result:.5g} V")

def proc_resistance_electricity():
    from physics.dimensional_formulas import resistance
    voltage = float(input("Enter voltage (V): "))
    current = float(input("Enter current (A): "))
    result = resistance(voltage, current)
    print(f"Resistance = {result:.5g} Ω")

def proc_capacitance():
    from physics.dimensional_formulas import capacitance
    q = float(input("Enter charge (C): "))
    voltage = float(input("Enter voltage (V): "))
    result = capacitance(q, voltage)
    print(f"Capacitance = {result:.5g} F")

def proc_conductance():
    from physics.dimensional_formulas import conductance
    r = float(input("Enter resistance (Ω): "))
    result = conductance(r)
    print(f"Conductance = {result:.5g} S")

def proc_electric_field():
    from physics.dimensional_formulas import electric_field
    F = float(input("Enter electric force (N): "))
    q = float(input("Enter charge (C): "))
    result = electric_field(F, q)
    print(f"Electric Field = {result:.5g} N/C")

def proc_dielectric_constant():
    from physics.dimensional_formulas import dielectric_constant
    eps_mat = float(input("Enter permittivity of material (F/m): "))
    eps0 = float(input("Enter permittivity of free space (F/m): "))
    result = dielectric_constant(eps_mat, eps0)
    print(f"Dielectric Constant = {result:.5g}")

def proc_specific_charge():
    from physics.dimensional_formulas import specific_charge
    q = float(input("Enter charge (C): "))
    m = float(input("Enter mass (kg): "))
    result = specific_charge(q, m)
    print(f"Specific Charge = {result:.5g} C/kg")


# ==========================================================
# MAGNETISM
# ==========================================================

def proc_magnetic_flux():
    from physics.dimensional_formulas import magnetic_flux
    B = float(input("Enter magnetic field (T): "))
    A = float(input("Enter area (m²): "))
    result = magnetic_flux(B, A)
    print(f"Magnetic Flux = {result:.5g} Wb")

def proc_inductance():
    from physics.dimensional_formulas import inductance
    phi = float(input("Enter magnetic flux (Wb): "))
    I = float(input("Enter current (A): "))
    result = inductance(phi, I)
    print(f"Inductance = {result:.5g} H")

def proc_magnetic_dipole_moment():
    from physics.dimensional_formulas import magnetic_dipole_moment
    tau = float(input("Enter torque (N·m): "))
    B = float(input("Enter magnetic field (T): "))
    result = magnetic_dipole_moment(tau, B)
    print(f"Magnetic Dipole Moment = {result:.5g} A·m²")

def proc_magnetic_field_strength():
    from physics.dimensional_formulas import magnetic_field_strength
    mu = float(input("Enter magnetic moment (A·m²): "))
    V = float(input("Enter volume (m³): "))
    result = magnetic_field_strength(mu, V)
    print(f"Magnetic Field Strength = {result:.5g} A/m")

def proc_magnetic_moment():
    from physics.dimensional_formulas import magnetic_moment
    I = float(input("Enter current (A): "))
    A = float(input("Enter area (m²): "))
    result = magnetic_moment(I, A)
    print(f"Magnetic Moment = {result:.5g} A·m²")

def proc_self_inductance_coefficient():
    from physics.dimensional_formulas import coefficient_self_induction
    emf = float(input("Enter induced emf (V): "))
    di_dt = float(input("Enter rate of change of current (A/s): "))
    result = coefficient_self_induction(emf, di_dt)
    print(f"Self-Induction Coefficient = {result:.5g} H")

def proc_mutual_inductance_coefficient():
    from physics.dimensional_formulas import coefficient_mutual_induction
    emf = float(input("Enter mutual emf (V): "))
    di_dt = float(input("Enter rate of change of current (A/s): "))
    result = coefficient_mutual_induction(emf, di_dt)
    print(f"Mutual Inductance Coefficient = {result:.5g} H")

def proc_relative_permeability():
    from physics.dimensional_formulas import relative_permeability
    mu_r = float(input("Enter permeability of material (H/m): "))
    mu0 = float(input("Enter permeability of free space (H/m): "))
    result = relative_permeability(mu_r, mu0)
    print(f"Relative Permeability = {result:.5g}")

def proc_magnetic_flux_density():
    from physics.dimensional_formulas import magnetic_flux_density
    I = float(input("Enter current (A): "))
    dl = float(input("Enter length element (m): "))
    r = float(input("Enter distance (m): "))
    result = magnetic_flux_density(I, dl, r)
    print(f"Magnetic Flux Density = {result:.5g} T")


__all__ = [name for name in globals() if not name.startswith("_")]