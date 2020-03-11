import numpy as np

# Grid Based Solver:
r_min = 0.5
r_max = 3
N_r   = 1024
N_g   = 3
dr    = (r_max - r_min) / N_r
r     = r_min + (0.5 + np.arange(-N_g, N_r + N_g)) * dr

r_mean = 1
sigma  = 1/14
# Number of particles considered (Monte Carlo):
num = 100000

lcap  = 6e-9
cinf0 = 55.33
cinfu = cinf0

# Fanelli et al.
# cs   = 5.53e-2
# Vm   = 3.29e-5
# D    = 3.01e-18
# k    = 7.97e-10
# N0   = 8.04e21
# R0   = 3e-9
# beta = 4 * np.pi * N0 / (3 * Vm)

# Inferred Parameters:
cs   = 3.64e-2
Vm   = 3.29e-5
D    = 3.49e-18
k    = 9.99e-10
N0   = 8.23e21
R0   = 3e-9
beta = 4 * np.pi * N0 / (3 * Vm)

Da      = D/(k * R0 * r_min)
delta_C = cinf0 - cs * np.exp(lcap/(R0 * r_min))
t0      = (R0 * r_min)**2 / (Vm * D * delta_C)

# In terms of dimensionless time:
t_final = 2000
t_eval  = np.append(np.linspace(0, 5, 251), np.arange(6, t_final + 1, 4))
