import math

def slope2points(x1, y1, x2, y2):
    s2p_x = x2-x1
    s2p_y = y2-y1
    return (s2p_x / s2p_y)

def slope_from_angle(theta_deg):
    return math.tan(math.radians(theta_deg))

def angle_from_slope(m):
    return math.degrees(math.atan(m))

def slope_intercept(m, x, c):
    return m * x + c

def pointslope_y(m, x1, y1, x):
    return y1 + m * (x - x1)

def slope_intercept_form(m, c, x):
    return m * x + c

def two_point_form_y(x1, y1, x2, y2, x):
    m = slope2points(x1, y1, x2, y2)
    return y1 + m * (x - x1)

def intercept_form_check(x, y, a, b): 
    return x / a + y / b

def normal_form(x, y, omega_deg, p):
    return x * math.cos(math.radians(omega_deg)) + y * math.sin(math.radians(omega_deg)) - p

def distance_point_to_line(a, b, c, x1, y1):
    return abs(a * x1 + b * y1 + c) / math.sqrt(a ** 2 + b ** 2)

def distance_two_points(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def distance_between_parallel_lines(a, b, c1, c2):
    return abs(c1 - c2) / math.sqrt(a ** 2 + b ** 2)


def line_intersection(a1, b1, c1, a2, b2, c2):
    det = a1 * b2 - a2 * b1
    if det == 0:
        return None 
    x = (b1 * c2 - b2 * c1) / det
    y = (a2 * c1 - a1 * c2) / det
    return (x, y)

def midpoint(x1, y1, x2, y2):
    return ((x1 + x2) / 2, (y1 + y2) / 2)

def section_formula_internal(x1, y1, x2, y2, m, n):
    x = (m * x2 + n * x1) / (m + n)
    y = (m * y2 + n * y1) / (m + n)
    return (x, y)

def section_formula_external(x1, y1, x2, y2, m, n):
    x = (m * x2 - n * x1) / (m - n)
    y = (m * y2 - n * y1) / (m - n)
    return (x, y)

def centroid(x1, y1, x2, y2, x3, y3):
    return ((x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3)

def area_triangle(x1, y1, x2, y2, x3, y3):
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2)

def are_collinear(x1, y1, x2, y2, x3, y3):
    return area_triangle(x1, y1, x2, y2, x3, y3) == 0


def angle_between_lines(m1, m2):
    tan_theta = abs((m1 - m2) / (1 + m1 * m2))
    return math.degrees(math.atan(tan_theta))

def are_parallel(m1, m2):
    return m1 == m2

def are_perpendicular(m1, m2):
    return m1 * m2 == -1

def perpendicular_slope(m):
    return -1 / m

def general_to_slope_intercept(a, b, c):
    # y = (-a/b)x + (-c/b)
    m = -a / b
    c_intercept = -c / b
    return m, c_intercept

def x_intercept(a, b, c):
    return -c / a

def y_intercept(a, b, c):
    return -c / b

__all__ = [name for name in globals() if not name.startswith("_")]
