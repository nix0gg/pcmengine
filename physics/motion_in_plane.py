import math
from physics_constants import g
def horizontal_displacement(u,theta,t):
    return u*math.cos(math.radians(theta))*t

def vertical_displacement(u,theta,t):
    return u*math.sin(math.radians(theta))*t - 0.5*g*t**2

def time_of_flight(u,theta):
    return (2*u*math.sin(math.radians(theta)))/g

def max_height(u,theta):
    return (u**2 * (math.sin(math.radians(theta)))**2) / (2*g)

def range_of_projectile(u,theta):
    return (u**2 * math.sin(math.radians(2*theta))) / g

def horizontal_velocity(u,theta):
    return u*math.cos(math.radians(theta))

def vertical_velocity(u,theta):
    return u*math.sin(math.radians(theta))

def resultant_velocity(u,theta,t):
    vx = horizontal_velocity(u,theta)
    vy = vertical_velocity(u,theta) - g*t
    return math.sqrt(vx**2 + vy**2)

def angle_of_velocity(u,theta,t):
    vx = horizontal_velocity(u,theta)
    vy = vertical_velocity(u,theta) - g*t
    return math.degrees(math.atan2(vy,vx))

__all__ = [name for name in globals() if not name.startswith("_")]
