import math
from physics_constants import g

def velocity(u, t, isDirectionUp=True):
    return u - g * t if isDirectionUp else u + g * t

def displacement(u, t, isDirectionUp=True):
    return u * t - 0.5 * g * t**2 if isDirectionUp else u * t + 0.5 * g * t**2

def velocity_from_displacement(u, h, isDirectionUp=True):
    v_squared = u**2 - 2 * g * h if isDirectionUp else u**2 + 2 * g * h
    return math.sqrt(v_squared) if v_squared >= 0 else float('nan')

def max_height(u):
    return u**2 / (2 * g)

def time_to_max_height(u):
    return u / g

def time_of_flight(u):
    return 2 * u / g

def ground_strike_velocity(h):
    return math.sqrt(2 * g * h)

def height_from_time(t):
    return 0.5 * g * t**2

def velocity_from_height(h):
    return math.sqrt(2 * g * h)

def time_from_height(h):
    return math.sqrt(2 * h / g)

def freefall_velocity(t):
    return g * t

def time_from_velocity(v):
    return v / g

def initial_velocity_from_displacement_final_velocity(v, h, isDirectionUp=True):
    u_squared = v**2 + 2 * g * h if isDirectionUp else v**2 - 2 * g * h
    return math.sqrt(u_squared) if u_squared >= 0 else float('nan')

def average_velocity(u, v):
    return (u + v) / 2

def average_acceleration(u, v, t):
    return (v - u) / t

def impulse_2(m, v, u):
    return m * (v - u)

def n1(m):
    return m * g

def n2(m, a):
    return m * (g + a)

def n3(m, a):
    return m * (g - a)

def m(n):
    return n / g

__all__ = [name for name in globals() if not name.startswith("_")]
