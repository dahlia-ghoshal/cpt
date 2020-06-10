import numpy as np

# physical constants (SI units)
e_charge = 1.602e-19
e_mass = 9.109e-31
epsilon_0 = 8.854e-12
hbar = 1.055e-34
bohr_rad = 4 * np.pi * epsilon_0 * hbar**2 / (e_mass * e_charge**2)
c = 2.998e8

# 87Rb constants
g_decay = 36.1  # [MHz] 87Rb D1
g31 = g_decay/2
g32 = g_decay/2
gf = g_decay/2  # decoherence rate between ground and excited states
gs = 0  # decoherence rate between ground states

n = 1  #refractive index
dipole = 2.992 * e_charge * bohr_rad  #transition dipole moment


def w1(I1):
    E1_amplitude = (2 / (n*epsilon_0*c) * I1) ** (1/2)
    w1 = dipole / hbar * E1_amplitude * 1e-6 # [MHz] Rabi frequency (assume light polarization parallel to dipole moment)
    return w1

def w2(I2):
    E2_amplitude = (2 / (n*epsilon_0*c) * I2) ** (1/2)
    w2 = dipole / hbar * E2_amplitude * 1e-6  # [MHz]
    return w2
    
    