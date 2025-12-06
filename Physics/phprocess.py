import os
import sys

try:
    import Physics.dimensional_formulas
    from Physics.physics_constants import g
except ModuleNotFoundError:
    # If running the file directly (not as a package), make project root importable
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    import Physics.dimensional_formulas
    from Physics.physics_constants import g

# ---------------------------- Basic Quantities ----------------------------
def proc_area():
    length = float(input("Enter length (L): "))
    breadth = float(input("Enter breadth (L): "))
    print("Area (L^2):", Physics.dimensional_formulas.area(length, breadth))

def proc_volume():
    length = float(input("Enter length (L): "))
    breadth = float(input("Enter breadth (L): "))
    height = float(input("Enter height (L): "))
    print("Volume (L^3):", Physics.dimensional_formulas.volume(length, breadth, height))

def proc_mass_density():
    mass = float(input("Enter mass (M): "))
    volume = float(input("Enter volume (L^3): "))
    print("Mass Density (M/L^3):", Physics.dimensional_formulas.mass_density(mass, volume))

def proc_frequency():
    time_period = float(input("Enter time period (T): "))
    print("Frequency (1/T):", Physics.dimensional_formulas.frequency(time_period))

def proc_speed():
    distance = float(input("Enter distance (L): "))
    time = float(input("Enter time (T): "))
    print("Speed (L/T):", Physics.dimensional_formulas.speed(distance, time))

def proc_acceleration():
    change_in_velocity = float(input("Enter change in velocity (L/T): "))
    time = float(input("Enter time (T): "))
    print("Acceleration (L/T^2):", Physics.dimensional_formulas.acceleration(change_in_velocity, time))

# ---------------------------- Mechanics ----------------------------
def proc_force():
    mass = float(input("Enter mass (M): "))
    acceleration = float(input("Enter acceleration (L/T^2): "))
    print("Force (MLT^-2):", Physics.dimensional_formulas.force(mass, acceleration))

def proc_impulse():
    force = float(input("Enter force (MLT^-2): "))
    time = float(input("Enter time (T): "))
    print("Impulse (MLT^-1):", Physics.dimensional_formulas.impulse(force, time))

def proc_work():
    force = float(input("Enter force (MLT^-2): "))
    distance = float(input("Enter distance (L): "))
    print("Work (ML^2T^-2):", Physics.dimensional_formulas.work(force, distance))

def proc_power():
    work = float(input("Enter work done (ML^2T^-2): "))
    time = float(input("Enter time taken (T): "))
    print("Power (ML^2T^-3):", Physics.dimensional_formulas.power(work, time))

def proc_momentum():
    mass = float(input("Enter mass (M): "))
    velocity = float(input("Enter velocity (L/T): "))
    print("Momentum (MLT^-1):", Physics.dimensional_formulas.momentum(mass, velocity))

def proc_pressure():
    force = float(input("Enter force (MLT^-2): "))
    area = float(input("Enter area (L^2): "))
    print("Pressure (ML^-1T^-2):", Physics.dimensional_formulas.pressure(force, area))

def proc_strain():
    change = float(input("Enter change in dimension (L): "))
    original = float(input("Enter original dimension (L): "))
    print("Strain (dimensionless):", Physics.dimensional_formulas.strain(change, original))

def proc_modulus_of_elasticity():
    stress = float(input("Enter stress (ML^-1T^-2): "))
    strain = float(input("Enter strain (dimensionless): "))
    print("Modulus of Elasticity (ML^-1T^-2):", Physics.dimensional_formulas.modulus_of_elasticity(stress, strain))

def proc_surface_tension():
    force = float(input("Enter force (MLT^-2): "))
    length = float(input("Enter length (L): "))
    print("Surface Tension (MT^-2):", Physics.dimensional_formulas.surface_tension(force, length))

def proc_kinetic_energy():
    mass = float(input("Enter mass (M): "))
    velocity = float(input("Enter velocity (L/T): "))
    print("Kinetic Energy (ML^2T^-2):", Physics.dimensional_formulas.kinetic_energy(mass, velocity))

def proc_potential_energy():
    mass = float(input("Enter mass (M): "))
    height = float(input("Enter height (L): "))
    print("Potential Energy (ML^2T^-2):", Physics.dimensional_formulas.potential_energy(mass, height, g))

def proc_escape_velocity():
    mass = float(input("Enter mass of body (M): "))
    radius = float(input("Enter radius (L): "))
    G = float(input("Enter gravitational constant (L^3M^-1T^-2): "))
    print("Escape Velocity (L/T):", Physics.dimensional_formulas.escape_velocity(mass, radius, G))

# ---------------------------- Electricity ----------------------------
def proc_charge_dimensional_formula():
    current = float(input("Enter current (A): "))
    time = float(input("Enter time (T): "))
    print("Charge (AT):", Physics.dimensional_formulas.charge(current, time))

def proc_current_density():
    current = float(input("Enter current (A): "))
    area = float(input("Enter area (L^2): "))
    print("Current Density (AL^-2):", Physics.dimensional_formulas.current_density(current, area))

def proc_voltage():
    work = float(input("Enter work done (ML^2T^-2): "))
    charge = float(input("Enter charge (AT): "))
    print("Voltage (ML^2T^-3A^-1):", Physics.dimensional_formulas.voltage(work, charge))

def proc_resistance():
    voltage = float(input("Enter voltage (ML^2T^-3A^-1): "))
    current = float(input("Enter current (A): "))
    print("Resistance (ML^2T^-3A^-2):", Physics.dimensional_formulas.resistance(voltage, current))

# ---------------------------- Optics and Radiation ----------------------------
def proc_refractive_index():
    speed_in_medium = float(input("Enter speed in medium (L/T): "))
    print("Refractive Index (dimensionless):", Physics.dimensional_formulas.refractive_index(speed_in_medium))

def proc_wavenumber():
    wavelength = float(input("Enter wavelength (L): "))
    print("Wavenumber (L^-1):", Physics.dimensional_formulas.wavenumber(wavelength))

# ---------------------------- Heat and Thermodynamics ----------------------------
def proc_heat_capacity():
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    temp_change = float(input("Enter temperature change (K): "))
    print("Heat Capacity (ML^2T^-2K^-1):", Physics.dimensional_formulas.heat_capacity(heat_energy, temp_change))

def proc_specific_heat_capacity():
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    mass = float(input("Enter mass (M): "))
    temp_change = float(input("Enter temperature change (K): "))
    print("Specific Heat Capacity (L^2T^-2K^-1):", 
          Physics.dimensional_formulas.specific_heat_capacity(heat_energy, mass, temp_change))

def proc_latent_heat():
    heat_energy = float(input("Enter heat energy (ML^2T^-2): "))
    mass = float(input("Enter mass (M): "))
    print("Latent Heat (L^2T^-2):", Physics.dimensional_formulas.latent_heat(heat_energy, mass))

# ---------------------------- Motion in One Dimension ----------------------------


def proc_velocity_1d():
    from Physics.motion_in_one_dimension import velocity
    u = float(input("Enter initial velocity (L/T): "))
    t = float(input("Enter time (T): "))
    direction = input("Is direction upward? (yes/no): ").lower()
    isDirectionUp = direction == "yes"
    print("Final Velocity (L/T):", velocity(u, t, isDirectionUp))

def proc_displacement_1d():
    from Physics.motion_in_one_dimension import displacement
    u = float(input("Enter initial velocity (L/T): "))
    t = float(input("Enter time (T): "))
    direction = input("Is direction upward? (yes/no): ").lower()
    isDirectionUp = direction == "yes"
    print("Displacement (L):", displacement(u, t, isDirectionUp))
    

def proc_max_height():
    from Physics.motion_in_one_dimension import max_height
    u = float(input("Enter initial velocity (L/T): "))
    print("Maximum Height (L):", max_height(u))

def proc_time_of_flight_1d():
    from Physics.motion_in_one_dimension import time_of_flight
    u = float(input("Enter initial velocity (L/T): "))
    print("Time of Flight (T):", time_of_flight(u))

def proc_ground_strike_velocity():
    from Physics.motion_in_one_dimension import ground_strike_velocity
    h = float(input("Enter height (L): "))
    print("Ground Strike Velocity (L/T):", ground_strike_velocity(h))

def proc_average_velocity():
    from Physics.motion_in_one_dimension import average_velocity
    u = float(input("Enter initial velocity (L/T): "))
    v = float(input("Enter final velocity (L/T): "))
    print("Average Velocity (L/T):", average_velocity(u, v))

def proc_average_acceleration():
    from Physics.motion_in_one_dimension import average_acceleration
    u = float(input("Enter initial velocity (L/T): "))
    v = float(input("Enter final velocity (L/T): "))
    t = float(input("Enter time (T): "))
    print("Average Acceleration (L/T^2):", average_acceleration(u, v, t))

def proc_normal_force():
    from Physics.motion_in_one_dimension import n1
    m = float(input("Enter mass (M): "))
    print("Normal Force (MLT^-2):", n1(m))

def proc_normal_force_accelerated():
    from Physics.motion_in_one_dimension import n2,n3
    m = float(input("Enter mass (M): "))
    a = float(input("Enter acceleration (L/T^2): "))
    up = input("Is acceleration upward? (yes/no): ").lower()
    if up == "yes":
        print("Normal Force (MLT^-2):", n2(m, a))
    else:
        print("Normal Force (MLT^-2):", n3(m, a))

# ---------------------------- Motion in Plane ----------------------------
def proc_projectile_range():
    from Physics.motion_in_plane import range_of_projectile
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Range of Projectile (L):", range_of_projectile(u, theta))

def proc_projectile_max_height():
    from Physics.motion_in_plane import max_height
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Maximum Height of Projectile (L):", max_height(u, theta))

def proc_projectile_time_of_flight():
    from Physics.motion_in_plane import time_of_flight
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    print("Time of Flight (T):", time_of_flight(u, theta))

def proc_horizontal_displacement():
    from Physics.motion_in_plane import horizontal_displacement
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Horizontal Displacement (L):", horizontal_displacement(u, theta, t))

def proc_vertical_displacement():
    from Physics.motion_in_plane import vertical_displacement
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Vertical Displacement (L):", vertical_displacement(u, theta, t))

def proc_resultant_velocity():
    from Physics.motion_in_plane import resultant_velocity
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Resultant Velocity (L/T):", resultant_velocity(u, theta, t))

def proc_angle_of_velocity():
    from Physics.motion_in_plane import angle_of_velocity
    u = float(input("Enter initial velocity (L/T): "))
    theta = float(input("Enter initial angle with horizontal (degrees): "))
    t = float(input("Enter time (T): "))
    print("Angle of Velocity (degrees):", angle_of_velocity(u, theta, t))

# ---------------------------- Work, Power and Energy ----------------------------


def proc_work_angular():
    from Physics.work_power_energy import work
    force = float(input("Enter force (MLT^-2): "))
    distance = float(input("Enter distance (L): "))
    theta = float(input("Enter angle between force and displacement (radians): "))
    print("Work Done (ML^2T^-2):", work(force, distance, theta))

def proc_total_mechanical_energy():
    from Physics.work_power_energy import total_mechanical_energy
    ke = float(input("Enter kinetic energy (ML^2T^-2): "))
    pe = float(input("Enter potential energy (ML^2T^-2): "))
    print("Total Mechanical Energy (ML^2T^-2):", total_mechanical_energy(ke, pe))

def proc_power_from_force():
    from Physics.work_power_energy import power_2
    force = float(input("Enter force (MLT^-2): "))
    velocity = float(input("Enter velocity (L/T): "))
    theta = float(input("Enter angle between force and velocity (radians): "))
    print("Power (ML^2T^-3):", power_2(force, velocity,theta))

def proc_work_energy_theorem():
    from Physics.work_power_energy import work_energy_theorem
    k2 = float(input("Enter final kinetic energy (ML^2T^-2): "))
    k1 = float(input("Enter initial kinetic energy (ML^2T^-2): "))
    print("Change in Kinetic Energy (ML^2T^-2):", work_energy_theorem(k2, k1))

def proc_efficiency():
    from Physics.work_power_energy import efficiency
    work_out = float(input("Enter work output (ML^2T^-2): "))
    work_in = float(input("Enter work input (ML^2T^-2): "))
    print("Efficiency (%):", efficiency(work_out, work_in))

def proc_gravitational_potential_energy():
    from Physics.work_power_energy import gravitational_potential_energy
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r1 = float(input("Enter initial separation (L): "))
    r2 = float(input("Enter final separation (L): "))
    print("Gravitational Potential Energy (ML^2T^-2):", 
          gravitational_potential_energy(m1, r1, r2, m2))

def proc_spring_potential_energy():
    from Physics.work_power_energy import spring_potential_energy
    k = float(input("Enter spring constant (MT^-2): "))
    x = float(input("Enter displacement from equilibrium (L): "))
    print("Spring Potential Energy (ML^2T^-2):", spring_potential_energy(k, x))

def proc_spring_kinetic_energy():
    from Physics.work_power_energy import kinetic_energy_spring
    k = float(input("Enter spring constant (MT^-2): "))
    x = float(input("Enter current displacement (L): "))
    xmax = float(input("Enter maximum displacement (L): "))
    print("Spring Kinetic Energy (ML^2T^-2):", kinetic_energy_spring(k, x, xmax))

# ---------------------------- Mechanical Properties of Solids ----------------------------
def proc_longitudinal_stress():
    from Physics.mechanical_properties_solids import longitudinal_stress
    force = float(input("Enter force (MLT^-2): "))
    area = float(input("Enter cross-sectional area (L^2): "))
    print("Longitudinal Stress (ML^-1T^-2):", longitudinal_stress(force, area))

def proc_hooke_law_stress():
    from Physics.mechanical_properties_solids import hooke_law_stress
    E = float(input("Enter value of 'E':"))
    strain = float(input("Enter value of strain:"))
    return hooke_law_stress(E,strain)

def proc_shearing_stress1():
    from Physics.mechanical_properties_solids import shearing_stress1
    F= float(input("Enter force:"))
    A = float(input("Enter area:"))
    return shearing_stress1(F,A) 

def proc_shearing_stress2():
    from Physics.mechanical_properties_solids import shearing_stress2
    eta = float(input("Enter eta:"))
    theta = float(input("Enter theta:"))
    return shearing_stress2(eta,theta)

def proc_shearing_strain():
    from Physics.mechanical_properties_solids import shearing_strain
    dx = float(input("Enter dx:"))
    L = float(input("Enter l:"))
    return shearing_strain(dx,L)

def proc_strain_longitudinal():
    from Physics.mechanical_properties_solids import longitudinal_strain
    dL = float(input("Enter change in length (m): "))
    L = float(input("Enter original length (m): "))
    result = longitudinal_strain(dL, L)
    print(f"Longitudinal Strain = {result:.5g}")

def proc_shear_modulus():
    from Physics.mechanical_properties_solids import shear_modulus
    f = float(input("Enter Force (N): "))
    a = float(input("Enter Area (m²): "))
    dx = float(input("Enter displacement (m): "))
    l = float(input("Enter original length (m): "))
    result = shear_modulus(f, a, dx, l)
    print(f"Shear Modulus = {result:.5g} N/m²")

def proc_bulk_modulus():
    from Physics.mechanical_properties_solids import bulk_modulus_elasticity2
    P = float(input("Enter Pressure (Pa): "))
    dV = float(input("Enter change in Volume (m³): "))
    V = float(input("Enter original Volume (m³): "))
    result = bulk_modulus_elasticity2(P, dV, V)
    print(f"Bulk Modulus = {result:.5g} N/m²")

def proc_poisson_ratio():
    from Physics.mechanical_properties_solids import poisson_ratio
    lateral_strain = float(input("Enter lateral strain: "))
    longitudinal_strain = float(input("Enter longitudinal strain: "))
    result = poisson_ratio(lateral_strain, longitudinal_strain)
    print(f"Poisson's Ratio = {result:.5g}")

def proc_poisson_ratio2():
    from Physics.mechanical_properties_solids import poissonratio2
    dr = float(input("Enter value for dr:"))
    r = float(input("Enter value for r:"))
    dl = float(input("Enter value for dl:"))
    l = float(input("Enter value for l:"))
    return poissonratio2(dr,r,dl,l)

def proc_elastic_potential_energy2():
    from Physics.mechanical_properties_solids import elastic_potential_energy2
    F = float(input("Enter Force (N): "))
    A = float(input("Enter Area (m²): "))
    dL = float(input("Enter change in length (m): "))
    L = float(input("Enter original length (m): "))
    result = elastic_potential_energy2(F, A, dL, L)
    print(f"Elastic Potential Energy = {result:.5g} J")

def proc_elastic_potential_energy1():
    from Physics.mechanical_properties_solids import elastic_potential_energy1
    
# ---------------------------- Gravitation ----------------------------
def proc_gravitational_force():
    from Physics.gravitation import universal_gravitational_force
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r = float(input("Enter distance between centers (L): "))
    print("Gravitational Force (MLT^-2):", universal_gravitational_force(m1, m2, r))

def proc_gravitational_acceleration():
    from Physics.gravitation import gravitational_acceleration
    m = float(input("Enter mass of the body (M): "))
    r = float(input("Enter distance from center (L): "))
    print("Gravitational Acceleration (LT^-2):", gravitational_acceleration(m, r))

def proc_orbital_velocity():
    from Physics.gravitation import orbital_velocity
    m = float(input("Enter mass of central body (M): "))
    r = float(input("Enter orbital radius (L): "))
    print("Orbital Velocity (LT^-1):", orbital_velocity(m, r))

def proc_orbital_period():
    from Physics.gravitation import time_period
    m = float(input("Enter mass of central body (M): "))
    r = float(input("Enter orbital radius (L): "))
    print("Orbital Period (T):", time_period(m, r))

def proc_weight():
    from Physics.gravitation import weight
    m = float(input("Enter mass (M): "))
    print("Weight (MLT^-2):", weight(m))

def proc_gravitational_pe():
    from Physics.gravitation import gravitational_potential_energy
    m1 = float(input("Enter mass of first body (M): "))
    m2 = float(input("Enter mass of second body (M): "))
    r = float(input("Enter distance between centers (L): "))
    print("Gravitational Potential Energy (ML^2T^-2):", gravitational_potential_energy(m1, m2, r))

# ==========================================================
#  MECHANICAL PROPERTIES OF FLUIDS
# ==========================================================

def proc_viscosity():
    from Physics.mechanical_properties_fluids import viscosity
    F = float(input("Enter Force (N): "))
    A = float(input("Enter Area (m²): "))
    dv_dy = float(input("Enter velocity gradient (s⁻¹): "))
    result = viscosity(F, A, dv_dy)
    print(f"Viscosity = {result:.5g} Pa·s")

def proc_bernoulli():
    from Physics.mechanical_properties_fluids import bernoulli_theorem
    P = float(input("Enter Pressure (Pa): "))
    rho = float(input("Enter Density (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    h = float(input("Enter height (m): "))
    v = float(input("Enter velocity (m/s): "))
    result = bernoulli_theorem(P, rho, g, h, v)
    print(f"Total Bernoulli Constant = {result:.5g} J/m³")

def proc_torricelli():
    from Physics.mechanical_properties_fluids import torricelli_law
    g = float(input("Enter g (m/s²): "))
    h = float(input("Enter height (m): "))
    p = float(input("Enter pressure p (Pa): "))
    p0 = float(input("Enter reference pressure p₀ (Pa): "))
    rho = float(input("Enter density (kg/m³): "))
    result = torricelli_law(g, h, p, p0, rho)
    print(f"Velocity (Torricelli's Law) = {result:.5g} m/s")

def proc_excess_pressure():
    from Physics.mechanical_properties_fluids import excess_pressure1
    S = float(input("Enter Surface Tension (N/m): "))
    R = float(input("Enter Radius (m): "))
    result = excess_pressure1(S, R)
    print(f"Excess Pressure = {result:.5g} Pa")

def proc_capillary_rise():
    from Physics.mechanical_properties_fluids import capillary_rise
    S = float(input("Enter Surface Tension (N/m): "))
    rho = float(input("Enter Density (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    r = float(input("Enter Radius of tube (m): "))
    theta = float(input("Enter contact angle θ (radians): "))
    result = capillary_rise(S, rho, g, r, theta)
    print(f"Capillary Rise = {result:.5g} m")

def proc_terminal_velocity():
    from Physics.mechanical_properties_fluids import terminal_velocity
    r = float(input("Enter Radius of sphere (m): "))
    rho = float(input("Enter density of fluid (kg/m³): "))
    rho0 = float(input("Enter density of body (kg/m³): "))
    g = float(input("Enter g (m/s²): "))
    eta = float(input("Enter Viscosity (Pa·s): "))
    result = terminal_velocity(r, rho, rho0, g, eta)
    print(f"Terminal Velocity = {result:.5g} m/s")

def proc_hydraulic_lift():
    from Physics.mechanical_properties_fluids import hydraulic_lift
    a2 = float(input("Enter Area₂ (m²): "))
    a1 = float(input("Enter Area₁ (m²): "))
    f1 = float(input("Enter Force₁ (N): "))
    result = hydraulic_lift(a2, a1, f1)
    print(f"Force on larger piston = {result:.5g} N")

def proc_absolute_pressure():
    from Physics.mechanical_properties_fluids import absolute_pressure
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
    from Physics.thermal_properties_of_matter import temperature_conversion1
    f = float(input("Enter temperature in Fahrenheit (°F): "))
    result = temperature_conversion1(f)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_temperature_from_kelvin():
    from Physics.thermal_properties_of_matter import temperature_conversion2
    k = float(input("Enter temperature in Kelvin (K): "))
    result = temperature_conversion2(k)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_temperature_from_rankine():
    from Physics.thermal_properties_of_matter import temperature_conversion3
    R = float(input("Enter temperature in Rankine (°R): "))
    result = temperature_conversion3(R)
    print(f"Temperature in Celsius = {result:.5g} °C")

def proc_thermal_linear_expansion():
    from Physics.thermal_properties_of_matter import thermal_linear_expansion1
    l0 = float(input("Enter original length (m): "))
    beta = float(input("Enter coefficient of linear expansion (1/°C): "))
    delta_T = float(input("Enter change in temperature (°C): "))
    result = thermal_linear_expansion1(l0, beta, delta_T)
    print(f"Final length = {result:.5g} m")

def proc_thermal_volume_expansion():
    from Physics.thermal_properties_of_matter import thermal_volume_expansion
    V0 = float(input("Enter original volume (m³): "))
    gamma = float(input("Enter coefficient of volume expansion (1/°C): "))
    delta_T = float(input("Enter temperature change (°C): "))
    result = thermal_volume_expansion(V0, gamma, delta_T)
    print(f"Final volume = {result:.5g} m³")

def proc_latent_heat_fusion():
    from Physics.thermal_properties_of_matter import latent_heat_fusion
    m = float(input("Enter mass (kg): "))
    lf = float(input("Enter latent heat of fusion (J/kg): "))
    result = latent_heat_fusion(m, lf)
    print(f"Heat absorbed/released = {result:.5g} J")

def proc_specific_heat_capacity_tpm():
    from Physics.thermal_properties_of_matter import specific_heat_capacity
    Q = float(input("Enter heat absorbed (J): "))
    m = float(input("Enter mass (kg): "))
    delta_T = float(input("Enter temperature change (°C): "))
    result = specific_heat_capacity(Q, m, delta_T)
    print(f"Specific Heat Capacity = {result:.5g} J/kg·°C")

def proc_rate_of_heat_flow():
    from Physics.thermal_properties_of_matter import rate_of_heat_flow1
    k = float(input("Enter thermal conductivity (W/m·K): "))
    A = float(input("Enter Area (m²): "))
    delta_T = float(input("Enter temperature difference (°C): "))
    d = float(input("Enter thickness (m): "))
    result = rate_of_heat_flow1(k, A, delta_T, d)
    print(f"Rate of Heat Flow = {result:.5g} W")

def proc_newton_law_cooling():
    from Physics.thermal_properties_of_matter import newton_law_cooling
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
    from Physics.dimensional_formulas import charge
    current = float(input("Enter current (A): "))
    time = float(input("Enter time (s): "))
    result = charge(current, time)
    print(f"Charge = {result:.5g} C")

def proc_current_density_electricity():
    from Physics.dimensional_formulas import current_density
    current = float(input("Enter current (A): "))
    area = float(input("Enter area (m²): "))
    result = current_density(current, area)
    print(f"Current Density = {result:.5g} A/m²")

def proc_voltage_electricity():
    from Physics.dimensional_formulas import voltage
    work = float(input("Enter work done (J): "))
    q = float(input("Enter charge (C): "))
    result = voltage(work, q)
    print(f"Voltage = {result:.5g} V")

def proc_resistance_electricity():
    from Physics.dimensional_formulas import resistance
    voltage = float(input("Enter voltage (V): "))
    current = float(input("Enter current (A): "))
    result = resistance(voltage, current)
    print(f"Resistance = {result:.5g} Ω")

def proc_capacitance():
    from Physics.dimensional_formulas import capacitance
    q = float(input("Enter charge (C): "))
    voltage = float(input("Enter voltage (V): "))
    result = capacitance(q, voltage)
    print(f"Capacitance = {result:.5g} F")

def proc_conductance():
    from Physics.dimensional_formulas import conductance
    r = float(input("Enter resistance (Ω): "))
    result = conductance(r)
    print(f"Conductance = {result:.5g} S")

def proc_electric_field():
    from Physics.dimensional_formulas import electric_field
    F = float(input("Enter electric force (N): "))
    q = float(input("Enter charge (C): "))
    result = electric_field(F, q)
    print(f"Electric Field = {result:.5g} N/C")

def proc_dielectric_constant():
    from Physics.dimensional_formulas import dielectric_constant
    eps_mat = float(input("Enter permittivity of material (F/m): "))
    eps0 = float(input("Enter permittivity of free space (F/m): "))
    result = dielectric_constant(eps_mat, eps0)
    print(f"Dielectric Constant = {result:.5g}")

def proc_specific_charge():
    from Physics.dimensional_formulas import specific_charge
    q = float(input("Enter charge (C): "))
    m = float(input("Enter mass (kg): "))
    result = specific_charge(q, m)
    print(f"Specific Charge = {result:.5g} C/kg")


# ==========================================================
# MAGNETISM
# ==========================================================

def proc_magnetic_flux():
    from Physics.dimensional_formulas import magnetic_flux
    B = float(input("Enter magnetic field (T): "))
    A = float(input("Enter area (m²): "))
    result = magnetic_flux(B, A)
    print(f"Magnetic Flux = {result:.5g} Wb")

def proc_inductance():
    from Physics.dimensional_formulas import inductance
    phi = float(input("Enter magnetic flux (Wb): "))
    I = float(input("Enter current (A): "))
    result = inductance(phi, I)
    print(f"Inductance = {result:.5g} H")

def proc_magnetic_dipole_moment():
    from Physics.dimensional_formulas import magnetic_dipole_moment
    tau = float(input("Enter torque (N·m): "))
    B = float(input("Enter magnetic field (T): "))
    result = magnetic_dipole_moment(tau, B)
    print(f"Magnetic Dipole Moment = {result:.5g} A·m²")

def proc_magnetic_field_strength():
    from Physics.dimensional_formulas import magnetic_field_strength
    mu = float(input("Enter magnetic moment (A·m²): "))
    V = float(input("Enter volume (m³): "))
    result = magnetic_field_strength(mu, V)
    print(f"Magnetic Field Strength = {result:.5g} A/m")

def proc_magnetic_moment():
    from Physics.dimensional_formulas import magnetic_moment
    I = float(input("Enter current (A): "))
    A = float(input("Enter area (m²): "))
    result = magnetic_moment(I, A)
    print(f"Magnetic Moment = {result:.5g} A·m²")

def proc_self_inductance_coefficient():
    from Physics.dimensional_formulas import coefficient_self_induction
    emf = float(input("Enter induced emf (V): "))
    di_dt = float(input("Enter rate of change of current (A/s): "))
    result = coefficient_self_induction(emf, di_dt)
    print(f"Self-Induction Coefficient = {result:.5g} H")

def proc_mutual_inductance_coefficient():
    from Physics.dimensional_formulas import coefficient_mutual_induction
    emf = float(input("Enter mutual emf (V): "))
    di_dt = float(input("Enter rate of change of current (A/s): "))
    result = coefficient_mutual_induction(emf, di_dt)
    print(f"Mutual Inductance Coefficient = {result:.5g} H")

def proc_relative_permeability():
    from Physics.dimensional_formulas import relative_permeability
    mu_r = float(input("Enter permeability of material (H/m): "))
    mu0 = float(input("Enter permeability of free space (H/m): "))
    result = relative_permeability(mu_r, mu0)
    print(f"Relative Permeability = {result:.5g}")

def proc_magnetic_flux_density():
    from Physics.dimensional_formulas import magnetic_flux_density
    I = float(input("Enter current (A): "))
    dl = float(input("Enter length element (m): "))
    r = float(input("Enter distance (m): "))
    result = magnetic_flux_density(I, dl, r)
    print(f"Magnetic Flux Density = {result:.5g} T")


__all__ = [name for name in globals() if not name.startswith("_")]