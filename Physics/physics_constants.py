import math

#Speed of light in vacuuum
c = 299792458  # m/s

#Planck constant
h = 6.62607015e-34 

# Reduced Planck constant
hbar = 1.054571817e-34  # J·s

#Boltzmann constant
kB = 1.380649e-23  # J/K

#Newtonian constant of gravitation
G = 6.67430e-11  # m^3·kg^−

#Cosmological constant
Lambda = 1.1056e-52  # m^−2

#Stefan-Boltzmann constant
sigma = 5.670374419e-8  # W·m^−2·K^−4

#First radiation constant
c1 = 3.741771852e-16  # W·m^2

#Second radiation constant
c2 = 1.438776877e-2  # m·K

#Wien wavelength displacement law constant
b = 2.897771955e-3  # m·K

#Wien entropy displacement law constant
b_entropy = 3.00292e-3  # m·K

#Elementary charge
ec = 1.602176634e-19  # C

#Conductance quantum
G0 = 7.748091729e-5  # S

#Inverse conductance quantum
G01 = 12906.40372547  # Ω

#Von Klitzing constant
R_k = 25812.807451  # Ω

#Josephson constant
KJ = 483597.848416  # GHz/V

#Magnetic flux quantum
phi0 = 2.067833848e-15  # Wb

#Fine structure constant
alpha = 7.2973525693e-3  # dimensionless

#Inverse fine structure constant
alpha_inv = 137.035999084  # dimensionless

#Vacuum magnetic permeability
mu0 = 1.25663706212e-6  # N/A^

#characteristic impendance of vacuum
ZZ0 = 376.730313668  # Ω

#Vacuum electric permittivity
epsilon0 = 8.8541878128e-12  # F/m

#Electron mass
me = 9.1093837015e-31  # kg

#Muon mass
mmu = 1.883531627e-28  # kg

#Tau mass
mtau = 3.16754e-27  # kg

#Proton mass
mp = 1.67262192369e-27  # kg

#Neutron mass
mn = 1.67492749804e-27  # kg

#Proton to electron mass ratio
mp_me = 1836.15267343  # dimensionless

#W to Z mass ratio
mWmz = 0.88153  # dimensionless

#Sine square weak mixing angle
sin2_theta_W = 0.22290  # dimensionless

#Electron g-factor
ge = -2.00231930436256  # dimensionless

#Muon g-factor
gmu = -2.0023318418  # dimensionless

#Proton g factor
gp = 5.585694713  # dimensionless

#Quantum of circulation
kappa_0 = 3.6369475516e-4  # m^2/s

#Bohr magneton
muB = 9.2740100783e-24  # J/T

#Nuclear magneton
muN = 5.0507837461e-27  # J/T

#Classical electron radius
re = 2.8179403262e-15  # m

#Thompson cross section
sigmaT = 6.6524587321e-29  # m^2

#Bohr radius
a0 = 5.29177210903e-11  # m

#Rydberg constant
R_inf = 10973731.568160  # m^−1

#Rydberg unit of energy
R_y = 2.1798723611035e-18  # J

#Hartree energy
E_h = 4.3597447222071e-18  # J

#Fermi coupling constant
G_F = 1.1663787e-5  # GeV^−2

#Avagadro constant
N_A = 6.02214076e23  # mol^−1

#Molar gas constant
R = 8.314462618  # J·mol^−1·K^

#Faraday constant
F = 96485.33212  # C·mol^−1

#Molar planck constant
N_A_h = 3.9903127176e-10  # J·s

#Molar mass of carbon-12
M_u = 1e-3  # kg·mol^−1

#Atomic mass constant
m_u = 1.66053906660e-27  # kg

#Molar mass constant
M_u = 1e-3  # kg·mol^−1

#Molar volume of silicon
V_m_Si = 12.0588349e-6  # m^3·mol^−1

#Hyperfine transition frequency of cesium-133
Delta_nu_Cs = 9.192631770e9  # Hz

#Acceleration due to gravity
g = 9.80665  # m/s^2

# Permittivity of free space
epsilon_0 = 8.854187817e-12  # F/m

#Permeability of free space
mu_0 = 4 * math.pi * 1e-7  # N/A^2

#Mechanical equivalent of heat
J = 4.186  # J/cal

#Standard atmospheric pressure
p_atm= 1.01325e5  # Pa

#Absolute zero in Celsius
T_0_C = -273.15  # °C

#Electron volt
eV = 1.602176634e-19  # J

#Unified atomic mass unit
u = 1.66053906660e-27  # kg

#Electron rest energy
E_e_rest = 8.1871057769e-14  # J

#Energy equivalent of 1 atomic mass unit
uc2 = 931.5 #MeV

#volume of ideal gas 
V_m_ideal = 22.710947  # m^3·mol^−1 at STP

__all__ = [name for name in globals() if not name.startswith("_")]
