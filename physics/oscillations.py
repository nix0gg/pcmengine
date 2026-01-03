import math
def frequency(T):
    freq = 1/T
    return T
#Misc
def shm_x(A,omega,t,phi):
    x = A*math.sin(omega*t+phi)
    return x

def shm_y(A,omega,t,phi):
    y = A*math.cos(omega*t+phi)
    return y

def phase(omega,t,phi):
    phase = (omega*t)+phi
    return phase

#Angular SHM
def displacment_x_angular(A,omega,t):
    x = A*math.cos(omega*t)
    return x

def displacement_y_angular(A,omega,t):
    y= A*math.cos(omega*t)
    return y

def torque_angular(C,theta):
    t = -C *theta
    return t

def velocity_angular(C,I):
    omega = math.sqrt(C/I)
    return omega

def potential_energy1_angular(C,theta):
    U = (1/2)*C*(theta**2)
    return U

def potential_energy2_angular(I,omega,theta):
    U = (1/2)*I*(omega**2)*(theta**2)
    return U

def kinetic_energy_angular(I,omega):
    K = (1/2)*I*(omega**2)
    return K

def total_mechanical_energy_angular(I,omega,thetaknot):
    E = (1/2)*I*(omega**2)*thetaknot
    return E

def acceleration_angular(C,theta,I):
    alpha = -C*theta/I
#Simple pendulum

def time_period_short_length(l,g):
    T = 2*math.pi*math.sqrt(l/g)
    return T

def time_period_large_length(l,g,Re):
    T = 2*math.pi*math.sqrt(1/(g*((1/l)+(1/Re))))
    return T

#Oscillation of pendulum in different situations
def osccilation_in_liquid1(T_dash,T):
    row = T_dash /T
    return row

def oscillation_in_liquid2(g,geff):
    row = math.sqrt(g/geff)
    return row

def oscillation_in_liquid3(p,row2):
    row = math.sqrt(p/(p-row2))
    return row

def pendulum_in_lift(l,g,a):
    T = 2*math.pi*math.sqrt(l/(g+a))
    return T

#Other types of SHM
def utube1(l,g):
    T = 2*math.pi*math.sqrt(l/2*g)
    return T

def  utube2(h,g):
    T = 2*math.pi*math.sqrt(h/g)
    return T


def ball_of_radius_r_T(R,r,g):
    T = 2*math.pi*math.sqrt((R-r)/g)
    return T

def motion_cylindrical_piston(M,h,P,A):
    T = 2*math.pi*math.sqrt((M*h)/(P*A))
    return T

def torsional_pendulum(l,C):
    T= 2*math.pi*math.sqrt(l/C)
    return T

#Linear SHM
def displacement_linear_Y(A,omega,t,phi):
    Y = A*math.sin(omega*t+phi)
    return Y

def displacement_linear_X(A,omega,t,phi):
    X = A*math.sin(omega*t+phi)
    return X

def velocity_linear(omega,A,y):
    v = omega*math.sqrt(A**2-y**2)
    return v

def acceleration_linear(omega,y):
    a = (-omega**2)*y
    return a

def potential_energy_linear(m,omega,y):
    U_lin = 1/2*m*(omega**2)*(y**2)
    return U_lin

def kinetc_energy_lin(m,omega,A,y):
    K = 1/2*m*(omega**2)*((A**2)-(y**2))
    return K

def total_mechanical_energy_linear(m,omega,A):
    E = 1/2*m*(omega**2)*(A**2)
    return E

#SHM in spring block system
def time_period(m,k):
    T = 2*math.pi*math.sqrt(m/k)
    return T

def frequency_shm_sbs(k,m):
    f = 1/(2*math.pi)*math.sqrt(k/m)
    return f

def two_masses_spring_connection(mew,k):
    T = 2*math.pi*math.sqrt(mew,k)
    return T

def mew(m1,m2):
    mew = m1*m2/(m1+m2)
    return mew

#Oscillation of Spring Combination
def series_combination(m,k1,k2):
    T = 2*math.pi*math.sqrt((m*(k1+k2))/(k1*k2))
    return T

def parallel_combination(m,k1,k2):
    T = 2*math.pi*math.sqrt(m/(k1+k2))
    return T

def twospring_connection(m,k1,k2):
    T = 2*math.pi*math.sqrt(m/(k1+k2))

