import math

#Speed of transverse waves
def speed_of_transverse_waves_solid(eta,p):
    v = math.sqrt(eta/p)
    return v

def speed_of_transverse_wave_stretchedspring(T,mew):
    v= math.sqrt(T/mew)
    return v

#Speed of longitudinal waves

def slw_solid(B,eta,row):
    v = math.sqrt((B+(4/3)*eta)/row)
    return v

def slow_fluid(B,row):
    v = math.sqrt(B/row)

def newton_speed_of_sound1(B_iso,row):
    v = math.sqrt(B_iso/row)
    return v

def newton_of_speed_of_sound2(P,row):
    v= math.sqrt(P/row)
    return v

#Beats
def v_beat1(v1,v2):
    v_beat = v1-v2
    return v_beat

def v_beat2(v1,v2):
    v_beat = v2-v1
    return v_beat

def v2_1(v1,v_beat):
    v2 = v1+v_beat
    return v2

def v2_2(v1,v_beat):
    v2 = v1-v_beat
    return v2

#Stationary waves
def wave_formation(a,t,lmbda,x,T):
    y = 2*a*math.sin((2*math.pi*t)/lmbda)*math.cos((2*math.pi*x)/T)
    return y

def frequency_string1(n,v,L):
    v2 = (n*v)/(2*L)
    return v

def frequency_string2(n,L,T,mew):
    v = (n/2*L)*math.sqrt(T/mew)
    return v

#Organ Pipe
#Open Organ Pipe
def v1(v,L):
    v1 = v/(2*L)
    return v1

def vn_open(n,v,L):
    vn_open = (n*v)/(2*L)
    return vn_open

#Close Pipe
def v1(v,L):
    v1 = v/(4*L)
    return v1

def vn_close(n,v):
    vn_close = (2*n-1)*v
    return vn_close 

#Displacement relation for a progressive waves

def displacement_of_particle1(A,omega,t,k,x):
    y =A*math.sin((omega*t)-(k*x))
    return y

def displacement_of_particle2(A,t,T,x,lmbda):
    y = A*math.sin(2*math.pi*((t/T)-(x/lmbda)))
    return y

def displacement_of_particle3(A,lmbda,v,t,x):
    y = A*math.sin((2*math.pi/lmbda)*(v*t-x))
    return y

def amplitude1(v):
    A = 2*math.pi*v
    return A

def amplitude2(T):
    A = 2*math.pi*(1/T)
    return A

def phase_change_time(T,dT):
    dPhi = (2*math.pi/T)*dT
    return dPhi

def phase_change_position(dX,lmbda):
    dPhi = -(2*math.pi/lmbda)*dX
    return dPhi

