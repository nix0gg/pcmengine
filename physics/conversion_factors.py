#Length
one_meter = 0 
one_gallon = 0
one_km = 0.6215 #mi
one_mile = 1.609 #km
one_inch = 2.54 #cm
one_yd = 0
one_ly = 9.461e15 #m
one_angstrom = 0.1e-9 #m
one_feet = 0

#Area
one_kmsq = 0
one_msq = 0
one_inchsq = 6.4516 #cm2
one_feetsq = 9.29 * 1e-2 #m2
one_acre = 43560 #feet^2
one_milesq = 0
#Volume
one_meter3 = 1e6 #cm3
one_litre = 0
one_gallon = 0
one_inch3 = 16.39 #cm3
one_feet3 = 0

#Speed
one_kmph = 0
one_mph = 0

#Magnetic Field
one_g = 1e-4 #T
one_t = 0 

#Angle and angular speed
pi_rad = 180 #degrees
one_radian = 57.30 #degrees
one_degree = 1.745e-2 #radians
one_rev_per_min = 0.1047 #radians per second
one_radian_per_second = 9.549 #rev per min

#Mass
one_tonne = 0
one_kg = 0
one_u = 0 
one_slug = 14.59 #kg

#Density 
one_gcm = 0

#Force
one_newton =0
one_lbf = 4.4482 #N
one_kgf = 2.2046 #N

#Time
one_hour = 0
one_day = 0
one_year = 0

#Time 
one_hour = 0
one_day = 0
one_year = 0

#Pressure
one_pascal = 1 #Nms-2
one_atm = 0
one_torr = 0

#Energy
one_ftlbf = 0
one_btu = 0
one_kwh = 3.6e6 #J
one_cal = 4.186 #J
one_l_atm = 0
one_ev = 1.602e-19 #J
one_uc2 = 931.50e6 #eV
one_erg = 1e-7 #J
#Power
one_hp = 0
one_watt = 0
one_btu_per_min = 17.58 #W

#Thermal conductivity
one_wmk = 6.938 #btu in/(ft2 hr °F)
one_btu_in_hft2_f =0.1441 #W/mK

#User defined functions for multiple conversions
def meterconvert(convertTo):
    if convertTo == "yd":
        return 1.0936  # yards
    elif convertTo == "ft":
        return 3.281  # feet
    elif convertTo == "in":
        return 39.37  # inches
    return None

def yardconvert(convertTo):
    if convertTo == "ft":
        return 3.0  # feet
    elif convertTo == "cm":
        return 91.44  # centimeters
    return None

# Area Conversions
def kmsqconvert(convertTo):
    if convertTo == "misq":
        return 0.3861  # square miles
    elif convertTo == "acres":
        return 247.1  # acres
    return None

def msqconvert(convertTo):
    if convertTo == "cmsq":
        return 1e4  # square centimeters
    return None

def milesqconvert(convertTo):
    if convertTo == "acres":
        return 460.0  # acres
    elif convertTo == "kmsq":
        return 2.590  # square kilometers
    return None

# Volume Conversions
def litreconvert(convertTo):
    if convertTo == "cm3":
        return 1000.0  # cubic centimeters
    elif convertTo == "m3":
        return 1e-3  # cubic meters
    return None

def gallonconvert(convertTo):
    if convertTo == "L":
        return 3.786  # liters
    elif convertTo == "qt":
        return 4.0  # quarts
    elif convertTo == "pt":
        return 8.0  # pints
    elif convertTo == "oz":
        return 128.0  # fluid ounces
    elif convertTo == "in3":
        return 231.0  # cubic inches
    return None

def feet3convert(convertTo):
    if convertTo == "in3":
        return 1728.0  # cubic inches
    elif convertTo == "L":
        return 28.32  # liters
    elif convertTo == "cm3":
        return 2.832e4  # cubic centimeters
    return None

# Speed Conversions
def kmphconvert(convertTo):
    if convertTo == "ms":
        return 0.2778  # meters per second
    elif convertTo == "mph":
        return 0.6215  # miles per hour
    return None

# Magnetic Field Conversions
def teslaconvert(convertTo):
    if convertTo == "wbm2":
        return 1.0  # weber per square meter
    elif convertTo == "G":
        return 1e4  # gauss
    return None

# Mass Conversions
def tonneconvert(convertTo):
    if convertTo == "kg":
        return 1000.0  # kilograms
    elif convertTo == "Mg":
        return 1.0  # megagrams
    return None

def uconvert(convertTo):
    if convertTo == "kg":
        return 1.6606e-27  # kilograms
    elif convertTo == "MeV/c2":
        return 931.50  # MeV/c^2
    return None

# Density Conversions
def gcmconvert(convertTo):
    if convertTo == "kgm3":
        return 1000.0  # kilograms per cubic meter
    elif convertTo == "kgL":
        return 1.0  # kilograms per liter
    return None

# Force Conversions
def newtonconvert(convertTo):
    if convertTo == "lbf":
        return 0.2248  # pound-force
    elif convertTo == "dyn":
        return 1e5  # dynes
    return None

# Time Conversions
def hourconvert(convertTo):
    if convertTo == "min":
        return 60.0  # minutes
    elif convertTo == "ks":
        return 3.6  # kiloseconds
    return None

def dayconvert(convertTo):
    if convertTo == "h":
        return 24.0  # hours
    elif convertTo == "min":
        return 1440.0  # minutes
    elif convertTo == "ks":
        return 86.4  # kiloseconds
    return None

def yearconvert(convertTo):
    if convertTo == "d":
        return 365.24  # days
    elif convertTo == "Ms":
        return 31.56  # megaseconds
    return None

# Pressure Conversions
def atmconvert(convertTo):
    if convertTo == "kPa":
        return 101.325  # kilopascals
    elif convertTo == "mmHg":
        return 1.01325  # millimeters of mercury
    elif convertTo == "inHg":
        return 29.9  # inches of mercury
    elif convertTo == "ftH2O":
        return 33.8  # feet of water
    return None

def torrconvert(convertTo):
    if convertTo == "mmHg":
        return 1.0  # millimeters of mercury
    elif convertTo == "Pa":
        return 133.32  # pascals
    return None

# Energy Conversions
def ftlbfconvert(convertTo):
    if convertTo == "J":
        return 1.356  # joules
    elif convertTo == "kWh":
        return 1.286e-3  # kilowatt-hours
    return None

def latmconvert(convertTo):
    if convertTo == "J":
        return 101.325  # joules
    elif convertTo == "cal":
        return 24.217  # calories
    return None

def btuconvert(convertTo):
    if convertTo == "ftlb":
        return 778.0  # foot-pounds
    elif convertTo == "cal":
        return 252.0  # calories
    elif convertTo == "J":
        return 1054.35  # joules
    return None

# Power Conversions
def hpconvert(convertTo):
    if convertTo == "ftlbfs":
        return 550.0  # foot-pounds force per second
    elif convertTo == "W":
        return 745.7  # watts
    return None

def wattconvert(convertTo):
    if convertTo == "hp":
        return 1.341e-3  # horsepower
    elif convertTo == "ftlbfs":
        return 0.7376  # foot-pounds force per second
    return None

__all__ = [name for name in globals() if not name.startswith("_")]



