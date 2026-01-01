import math
from physics_constants import c, N_A, eV, ec

# ---------------------------- Basic Quantities ----------------------------

def area(length, breadth):
    return length * breadth

def volume(length, breadth, height):
    return length * breadth * height

def mass_density(mass, volume):
    return mass / volume

def frequency(time_period):
    return 1 / time_period

def speed(distance, time):
    return distance / time

def acceleration(change_in_velocity, time):
    return change_in_velocity / time


# ---------------------------- Mechanics ----------------------------

def force(mass, acceleration):
    return mass * acceleration

def impulse(force, time):
    return force * time

def work(force, distance):
    return force * distance

def power(work, time):
    return work / time

def momentum(mass, velocity):
    return mass * velocity  

def pressure(force, area):
    return force / area

def strain(change_in_dimension, original_dimension):
    return change_in_dimension / original_dimension

def modulus_of_elasticity(stress, strain):
    return stress / strain

def surface_tension(force, length):
    return force / length

def surface_energy(energy, area):
    return energy / area

def velocity_gradient(change_in_velocity, distance):
    return change_in_velocity / distance

def pressure_gradient(change_in_pressure, distance):
    return change_in_pressure / distance

def pressure_energy(pressure, volume):
    return pressure * volume

def coefficient_of_viscosity(shear_stress, velocity_gradient):
    return shear_stress / velocity_gradient

def angular_displacement(arc, radius):
    return arc / radius

def angular_velocity(angular_displacement, time):
    return angular_displacement / time

def angular_acceleration(angular_velocity, time):
    return angular_velocity / time

def moment_of_inertia(mass, radius_of_gyration):
    return mass * (radius_of_gyration ** 2)

def angular_momentum(moment_of_inertia, angular_velocity):
    return moment_of_inertia * angular_velocity

def torque(force, distance):
    return force * distance

def angular_frequency(time_period):
    return 2 * math.pi / time_period

def hubble_constant(recession_speed, distance):
    return recession_speed / distance

def intensity_of_wave(energy, time, area):
    return (energy / time) / area

def radiation_pressure(intensity, speed_of_light):
    return intensity / speed_of_light

def energy_density(energy, volume):
    return energy / volume

def critical_velocity(reynolds_no, viscosity_coefficient, density, radius):
    return (reynolds_no * viscosity_coefficient) / (density * radius)

def escape_velocity(mass, radius, gravitational_constant):
    return math.sqrt((2 * gravitational_constant * mass) / radius)

def kinetic_energy(mass, velocity):
    return 0.5 * mass * velocity ** 2

def potential_energy(mass, height, g):
    return mass * g * height

def rotational_kinetic_energy(moment_of_inertia, angular_velocity):
    return 0.5 * moment_of_inertia * angular_velocity ** 2

def efficiency(output_energy, input_energy):
    return (output_energy / input_energy) * 100

def angular_impulse(torque, time):
    return torque * time

def gravitational_constant(force, mass1, mass2, distance):
    return (force * distance ** 2) / (mass1 * mass2)

def gc2(f,r,m):
    return (f*r**2)/(m**2)

def poisson_ratio(lateral_strain, longitudinal_strain):
    return lateral_strain / longitudinal_strain

def couple(force,distance):
    return force * distance

def spring_constant(force, length):
    return force / length
                    
# ---------------------------- Heat and Thermodynamics ----------------------------

def plancks_constant(energy, frequency):
    return energy / frequency

def heat_capacity(heat_energy, temperature_change):
    return heat_energy / temperature_change

def specific_heat_capacity(heat_energy, mass, temperature_change):
    return heat_energy / (mass * temperature_change)

def latent_heat(heat_energy, mass):
    return heat_energy / mass

def thermal_expansion(change_in_dimension, original_dimension, temperature_change):
    return change_in_dimension / (original_dimension * temperature_change)

def bulk_modulus(volume, change_in_volume, change_in_pressure):
    return -change_in_pressure * volume / change_in_volume

def centripetal_acceleration(velocity, radius):
    return velocity ** 2 / radius

def stefan_constant(energy, temperature, area, time):
    return energy / (area * time * (temperature ** 4))

def wien_constant(wavelength, temperature):
    return wavelength * temperature

def boltzmann_constant(energy, temperature):
    return energy / temperature

def universal_gas_constant(pressure, volume, temperature, moles):
    return (pressure * volume) / (temperature * moles)

def mechanical_equivalence(w, h):
    return w / h

def specific_heat(heat_energy, mass, temperature_change):
    return heat_energy / (mass * temperature_change)

def thermal_capacity(m, specific_heat):
    return m * specific_heat

def molar_specific_heat(specific_heat, molar_mass):
    return specific_heat * molar_mass

def solar_constant(energy,time,area):
    return energy / (time * area)

def latent_heat(heat_energy, mass):
    return heat_energy / mass

def entropy(heat,temperature):
    return heat / temperature

def vanderwaals_constant_a(pressure, volume):
    return pressure * (volume ** 2)
def vanderwaals_constant_b(volume, b):
    b = volume
    return b

def thermal_conductivity(heat_energy, area, d0, dx, time):
    return heat_energy * (area * (d0 / dx) * time)

def temperature_gradient(d0, dx):
    return d0 / dx

def thermal_resistance(plane_thickness, thermal_cond, area):
    return plane_thickness / (thermal_cond * area)

def emissive_power(heat, area, time):
    return heat / (area * time)

def rate_of_cooling(lost_heat, time):
    return lost_heat / time

def molar_volume(volume, moles):
    return volume / moles
# ---------------------------- Electricity ----------------------------

def charge(current, time):
    return current * time

def current_density(current, area):
    return current / area

def voltage(work, charge):
    return work / charge

def resistance(voltage, current):
    return voltage / current

def capacitance(charge, voltage):
    return charge / voltage

def electrical_resistivity(resistance, length, area):
    return (resistance * area) / length

def electric_field(electric_force, charge):
    return electric_force / charge

def electric_flux(electric_field, area):
    return electric_field * area

def electric_dipole_moment(torque, electric_field):
    return torque / electric_field

def electric_intensity(potential_difference, distance):
    return potential_difference / distance
 
def electrical_permittivity(charge1, charge2, electric_force, distance):
    return (charge1 * charge2) / (4 * math.pi * electric_force * distance ** 2) #electrical permittivity = e0

def dielectric_constant(electrical_permittivity_material, electrical_permittivity_free_space):
    return electrical_permittivity_material / electrical_permittivity_free_space

def conductance(resistance):
    return 1 / resistance

def conductivity(resistivity):
    return 1 / resistivity

def electric_flux(electric_intensity, area):
    return electric_intensity * area

def specific_charge(charge, mass):
    return charge / mass
# ---------------------------- Magnetism ----------------------------

def magnetic_flux(magnetic_field, area):
    return magnetic_field * area

def inductance(magnetic_flux, current):
    return magnetic_flux / current

def magnetic_dipole_moment(torque, magnetic_field):
    return torque / magnetic_field

def magnetic_field_strength(magnetic_moment, volume):
    return magnetic_moment / volume

def permittivity_constant(charge1, charge2, electric_force, distance):
    return (charge1 * charge2) / (4 * math.pi * electric_force * distance ** 2)

def permeability_constant(force, distance, current1, current2, length):
    return (2 * math.pi * force * distance) / (current1 * current2 * length)

def magnetic_moment(current, area):
    return current * area

def pole_strength(magnetic_force, magnetic_field):
    return magnetic_force / magnetic_field

def coefficient_self_induction(induced_emf, rate_of_change_current):
    return induced_emf / rate_of_change_current

def coefficient_mutual_induction(mutual_emf, rate_of_change_current):
    return mutual_emf / rate_of_change_current

def magnetic_flux_density(I, dl, r):
    return (I*dl)/(4*math.pi*r**2)

def relative_permeability(permeability_material, permeability_free_space):
    return permeability_material / permeability_free_space

def magnetic_suspension_force(magnetic_field, pole_strength):
    return magnetic_field * pole_strength



# ---------------------------- Optics and Radiation ----------------------------

def refractive_index(speed_in_medium):
    return c / speed_in_medium

def faraday_constant():
    return N_A * ec

def wavenumber(wavelength):
    return (2 * math.pi) / wavelength

def radiant_flux(energy_emitted, time):
    return energy_emitted / time

def radiant_intensity(radiant_power, solid_angle):
    return radiant_power / solid_angle

def luminous_power(luminous_energy, time):
    return luminous_energy / time

def luminous_intensity(luminous_flux, solid_angle):
    return luminous_flux / solid_angle

def illumination_intensity(luminous_intensity, distance):
    return luminous_intensity / (distance ** 2)

def luminous_efficiency(total_luminous_flux, total_radiant_flux):
    return total_luminous_flux / total_radiant_flux

def illumination(luminous_flux_incident, area):
    return luminous_flux_incident / area

def luminous_flux(light_energy, seconds):
    return light_energy / seconds


# ---------------------------- Modern Physics ----------------------------

def mass_defect(sum_of_nucleon_masses, actual_nucleus_mass):
    return sum_of_nucleon_masses - actual_nucleus_mass

def binding_energy(mass_defect):
    return mass_defect * (c ** 2)


# ---------------------------- Oscillations & Circuits ----------------------------

def resonant_freq(inductance, capacitance):
    return 1 / (2 * math.pi * math.sqrt(inductance * capacitance))

def quality_factor_coil(resonant_freq, inductance, resistance):
    return (resonant_freq * inductance) / resistance

def power_lens(focal_length):
    return 1 / focal_length

def magnification(image_distance, object_distance):
    return image_distance / object_distance

def fluid_flow_rate(pressure, radius, viscosity_coefficient, length):
    return ((math.pi / 8) * pressure * (radius ** 4)) / (viscosity_coefficient * length)

def capacitive_reactance(angular_freq, capacitance):
    return 1 / (angular_freq * capacitance)

def inductive_reactance(angular_freq, inductance):
    return angular_freq * inductance

__all__ = [name for name in globals() if not name.startswith("_")]


