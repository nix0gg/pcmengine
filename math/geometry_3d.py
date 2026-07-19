import math

def distance(x1, y1, z1, x2, y2, z2):
    x= (x2-x1)**2
    y =(y2-y1)**2
    z= (z2-z1)**2 
    return math.sqrt(x+y+z)

def midpoint(x1, y1, z1, x2, y2, z2):
    midpoint_x=(x1+x2)/2
    midpoint_y=(y1+y2)/2
    midpoint_z=(z1+z2)/2
    return (midpoint_x,midpoint_y,midpoint_z)

def section_internal(x1, y1, z1, x2, y2, z2, m, n):
    x = (m*x2 + n*x1) / (m+n)
    y = (m*y2 + n*y1) / (m+n)
    z = (m*z2 + n*z1) / (m+n)
    return (x, y, z)

def section_external(x1, y1, z1, x2, y2, z2, m, n):
    x = (m*x2 - n*x1) / (m-n)
    y = (m*y2 - n*y1) / (m-n)
    z = (m*z2 - n*z1) / (m-n)
    return (x, y, z)

def direction_cosines(a, b, c):
    mag = math.sqrt(a**2 + b**2 + c**2)
    return (a/mag, b/mag, c/mag)

def direction_ratios_from_two_points(x1, y1, z1, x2, y2, z2):
    x = x2-x1
    y= y2-y1
    z=z2-z1
    return (x,y,z)

def dc_from_dr(a, b, c):
    return direction_cosines(a, b, c)

def sum_of_squares_dc(l, m, n):
    sum = (l**2)+(m**2)+(n**2)
    return sum

def angle_between_lines_dc(l1, m1, n1, l2, m2, n2):
    cos_theta = l1*l2 + m1*m2 + n1*n2
    return math.degrees(math.acos(cos_theta))

def angle_between_lines_dr(a1, b1, c1, a2, b2, c2):
    dot = a1*a2 + b1*b2 + c1*c2
    mag1 = math.sqrt(a1**2 + b1**2 + c1**2)
    mag2 = math.sqrt(a2**2 + b2**2 + c2**2)
    return math.degrees(math.acos(dot / (mag1 * mag2)))

def are_perpendicular_dr(a1, b1, c1, a2, b2, c2):
    a= a1*b2
    b= b1*b2
    c=c1*c2
    return(a+b+c==0)

def are_parallel_dr(a1, b1, c1, a2, b2, c2):
    a= a1/a2
    b= b1/b2
    c=c1/c2
    return(a==b==c)

def centroid_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3):
    c3d_x=(x1+x2+x3)/3
    c3d_y =(y1+y2+y3)/3
    c3d_z =  (z1+z2+z3)/3
    return (c3d_x, c3d_y ,c3d_z)

__all__ = [name for name in globals() if not name.startswith("_")]
