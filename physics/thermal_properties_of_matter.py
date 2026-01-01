def temperature_conversion1(f):
   return (f - 32) * 5.0 / 9.0

def temperature_conversion2(k):
    return k - 273.15

def temperature_conversion3(R):
     return (R - 491.67) * 5.0 / 9.0

def thermal_volume_expansion(V0,gamma,dk):
    return V0*(1+(gamma*dk))

def thermal_superficial_expansion(A0,beta,delta_T):
    return A0*(1+(beta*delta_T))

def thermal_linear_expansion1(l0,beta,delta_T):
    return l0*(1+beta*delta_T)

def thermal_linear_expansion2(beta):
    return beta / 2

def thermal_linear_expansion3(gamma):
    return gamma / 3

def latent_heat_fusion(m,lf):
    return m * lf

def latent_heat_vapourisation(m,lv):
    return m * lv

def latent_heat_sublimation(m,ls):
    return m * ls

def specific_heat_capacity(Q,m,delta_T):
    return Q / (m * delta_T)

def rate_of_heat_flow1(k,A,delta_T,d):
    return (k * A * delta_T) / d

def rate_of_flow_heat2(K,A,delta_T,dX):
    return K * A * (delta_T / dX)

def molar_heat_capacity(dQ,n,delta_T):
    return dQ / (n * delta_T)

def newton_law_cooling(k,m,s,t2,t1):
    return k * m * s * (t2 - t1)

def kirchhoff_law1(e1,a1):
    return e1 / a1

def kirchhoff_law2(e2,a2):
    return e2 / a2

def kirchhoff_law3(E, A):
    return E / A

__all__ = [name for name in globals() if not name.startswith("_")]
