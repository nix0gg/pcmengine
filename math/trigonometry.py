import math 

def degrees_to_radian(degrees):
    return degrees*(math.pi/180)

def radian_to_degree(radians):
    return radians *(180/math.pi)

def arc_length(r,theta_rad):
    return r*theta_rad

def area_sector(r,thetha_rad):
    return (0.5)*r**2*thetha_rad

def sin_d(degree):
    return math.sin(degrees_to_radian(degree))

def cos_d(degree):
    return math.cos(degrees_to_radian(degree))

def tan_d(degree):
    return math.tan(degrees_to_radian(degree))

def cosec_d(degree):
    return 1 / sin_d(degree)

def sec_d(degree):
    return 1 / cos_d(degree)

def cot_d(degree):
    return 1 / tan_d(degree)

def arcsin_d(x):
    return radian_to_degree(math.asin(x))

def arccos_d(x):
    return radian_to_degree(math.acos(x))

def arctan_x(x):
    return radian_to_degree(math.atan(x))

def sin_from_cos(val_cos):
    return math.sqrt(1-val_cos**2)

def cos_from_sin(val_sin):
    return math.sqrt(1-val_sin**2)

def sec_from_tan(val_tan):
    return math.sqrt(1+val_tan **2)

def cosec_from_cot(val_cot):
    return math.sqrt(1+val_cot**2)

def sin_sum(a,b):
    return sin_d(a)*cos_d(b) + cos_d(a)*sin_d(b)

def sin_diff(a,b):
    return cos_d(a)*cos_d(b) - sin_d(a)*sin_d(b)

def cos_sum(a,b):
    return cos_d(a)*cos_d(b)-sin_d(a)*sin_d(b)

def cos_diff(a,b):
    return cos_d(a)*cos_d(b)+sin_d(a)*sin_d(b)

def tan_sum(a,b):
    return (tan_d(a)+tan_d(b)) / (1-tan_d(a)**2)

def tan_diff(a,b):
    return ((tan_d(a)- tan_d(b))/1+tan_d(a)**2)

def sin_double(a):
    return 2*sin_d(a)*cos_d(a)

def cos_double1(a):
    return cos_d(a)**2-sin_d(a)**2

def cos_double2(a):
    return 1-2*sin_d(a)**2

def cos_double_3(a):
    return 2*cos_d(a)**2-1

def tan_double(a):
    return ((2*tan_d(a))/(1-tan_d(a)**2))

def sin_triple(a):
    return 3*sin_d(a) - 4*sin_d(a)**3

def cos_triple(a):
    return 4*cos_d(a)**3-3*cos_d(a)

def tan_triple(a):
    return (3*tan_d(a)-tan_d(a)**3)/(1-3*tan_d(a)**2)

def sin_half(a):
    return math.sqrt((1-cos_d(a))/2)

def cos_half(a):
    return math.sqrt((1+cos_d(a))/2)

def tan_half(a):
    return math.sqrt((1-cos_d(a))/2)

def sinxcos(a,b):
    return 0.5*(sin_sum(a,b)+sin_diff(a,b))

def cosxcos(a,b):
    return 0.5*(cos_diff(a,b)+cos_sum(a,b))

def sinxsin(a,b):
    return 0.5*(cos_diff(a,b) - cos_sum(a,b))

def sin_plus_sin(A, B): 
    return 2 * sin_d((A + B) / 2) * cos_d((A - B) / 2)

def sin_minus_sin(A, B): 
    return 2 * cos_d((A + B) / 2) * sin_d((A - B) / 2)

def cos_plus_cos(A, B): 
    return 2 * cos_d((A + B) / 2) * cos_d((A - B) / 2)

def cos_minus_cos(A, B): 
    return -2 * sin_d((A + B) / 2) * sin_d((A - B) / 2)

def sine_rule_side(a,A,B):
    return a*sin_d(B) / sin_d(A)

def cosine_rule_side(b,c,A):
    return math.sqrt(b**2+c**2-2*b*c*cos_d(A))

def cosine_rule_angle(a,b,c):
    return arccos_d((b**2+c**2-a**2)/(2*b*c))

def area_traingle_trig(b,c,A):
    return 0.5*b*c*sin_d(A)

__all__=[name for name in globals() if not name.startswith("_")]