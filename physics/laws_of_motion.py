import math
g = 9.8

def force(m, a):
    return m * a

def friction(mu, N):
    return mu * N

def inclined_plane_acceleration(theta):
    return g * math.sin(math.radians(theta))

def tension(m1, m2, theta):
    a = (m2 * g - m1 * g * math.sin(math.radians(theta))) / (m1 + m2)
    return m2 * (g - a)

def impulse(F, t):
    return F * t

def momentum(m, v):
    return m * v

def newtons_second_law(F, m):
    return F / m 

__all__ = [name for name in globals() if not name.startswith("_")]
