import math

def circle_equation(h, k, r, x, y):
    return (x - h) ** 2 + (y - k) ** 2 - r ** 2

def circle_radius(h, k, c):
    g, f = h, k
    return math.sqrt(g ** 2 + f ** 2 - c)

def circle_centre_general(g, f):
    return (-g, -f)

def circle_radius_general(g, f, c):
    return math.sqrt(g ** 2 + f ** 2 - c)

def circle_area(r):
    return math.pi * r ** 2

def circle_circumference(r):
    return 2 * math.pi * r

def parabola_focus_y2_4ax(a):
    return (a, 0)

def parabola_directrix_y2_4ax(a):
    return f"x = {-a}"

def parabola_latus_rectum(a):
    return 4 * a

def parabola_vertex():
    return (0, 0)

def parabola_y_from_x(a, x):
    return math.sqrt(4 * a * x)

def parabola_x2_4ay_focus(a):
    return (0, a)

def parabola_x2_4ay_directrix(a):
    return f"y = {-a}"


def ellipse_c(a, b):
    return math.sqrt(a ** 2 - b ** 2)

def ellipse_eccentricity(a, b):
    c = ellipse_c(a, b)
    return c / a

def ellipse_foci(a, b):
    c = ellipse_c(a, b)
    return ((-c, 0), (c, 0))

def ellipse_vertices(a):
    return ((-a, 0), (a, 0))

def ellipse_covertices(b):
    return ((0, -b), (0, b))

def ellipse_latus_rectum(a, b):
    return 2 * b ** 2 / a

def ellipse_directrix(a, b):
    e = ellipse_eccentricity(a, b)
    return a / e

def ellipse_area(a, b):
    return math.pi * a * b

def ellipse_semi_latus_rectum(a, b):
    return b ** 2 / a

def hyperbola_c(a, b):
    return math.sqrt(a ** 2 + b ** 2)

def hyperbola_eccentricity(a, b):
    c = hyperbola_c(a, b)
    return c / a

def hyperbola_foci(a, b):
    c = hyperbola_c(a, b)
    return ((-c, 0), (c, 0))

def hyperbola_vertices(a):
    return ((-a, 0), (a, 0))

def hyperbola_asymptotes(a, b):
    return (f"y = {b/a}x", f"y = {-b/a}x")

def hyperbola_latus_rectum(a, b):
    return 2 * b ** 2 / a

def hyperbola_directrix(a, b):
    e = hyperbola_eccentricity(a, b)
    return a / e


def conic_type(A, B, C):
    delta = B ** 2 - 4 * A * C
    if delta < 0:
        return "Ellipse (or circle if A=C, B=0)"
    elif delta == 0:
        return "Parabola"
    else:
        return "Hyperbola"

__all__ = [name for name in globals() if not name.startswith("_")]
