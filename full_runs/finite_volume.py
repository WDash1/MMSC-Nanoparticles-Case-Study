import numpy as np
from scipy.integrate import odeint
import h5py

lcap  = 6e-9
cinf0 = 55.33
cinfu = cinf0

cs   = 5.53e-2
Vm   = 3.29e-5
D    = 3.01e-18
k    = 7.97e-10
N0   = 8.04e21
R0   = 3e-9
beta = 4 * np.pi * N0 / (3 * Vm)
Da   = D/(k * R0)

delta_C = cinf0 - cs * np.exp(lcap/R0)
t0      = (R0**2 / (Vm * D * delta_C))

r_min = 0.1
r_max = 5
N_r   = 256
N_g   = 3
dr    = (r_max - r_min) / N_r

r  = r_min + (0.5 + np.arange(-N_g, N_r + N_g)) * dr
f0 = 10 * np.sqrt(0.5/np.pi) * np.exp(-50 * (r - 1)**2)

def minmod(x, y, z):

    min_of_all = np.minimum(np.minimum(abs(x), abs(y)), abs(z))

    signx = np.sign(x)
    signy = np.sign(y)
    signz = np.sign(z)
    
    result = 0.25 * abs(signx + signy) * (signx + signz) * min_of_all
    return result

def slope_minmod(input_array):
    """
    Reconstructs the input array using a 
    minmod limiter.

    Parameters
    ----------
    
    input_array: np.array
                 Array holding the cells data.
    """
    f_i_plus_one  = np.roll(input_array, -1)
    f_i_minus_one = np.roll(input_array,  1)
  
    forward_diff  = (f_i_plus_one - input_array  )
    backward_diff = (input_array  - f_i_minus_one)
    central_diff  = backward_diff + forward_diff

    # Sort of user controlled:
    slope_lim_theta = 2

    left   = slope_lim_theta * backward_diff
    center = 0.5 * central_diff
    right  = slope_lim_theta * forward_diff

    return(minmod(left, center, right))

def reconstruct_minmod(variable_to_reconstruct):
    
    slope = slope_minmod(variable_to_reconstruct)

    left_face_value  = variable_to_reconstruct - 0.5 * slope
    right_face_value = variable_to_reconstruct + 0.5 * slope
   
    return(left_face_value, right_face_value)

def riemann_solver(left_state, right_state, velocity):
    upwind_state = np.select([velocity >= 0, velocity < 0], 
                             [left_state, right_state]
                            )
    return upwind_state

def df_dt(f, t):
    cinfu   = cinf0 + beta * R0**3 * np.sum(f[:N_g] * r[:N_g]**3) * dr
    f[:N_g] = 0
    if(np.sum(f[-N_g:] * r[-N_g:]**3) * dr > 1e-8):
        raise Exception('Increase r_max!')
    
    f_left_plus_eps, f_right_minus_eps = reconstruct_minmod(f)
    f_left_minus_eps = np.roll(f_right_minus_eps, 1)

    cinf     = cinfu - beta * R0**3 * np.sum(f[N_g:-N_g] * r[N_g:-N_g]**3) * dr
    velocity = (cinf - cs * np.exp(lcap / (R0 * r))) / (delta_C * (Da + r))
    f_left   = riemann_solver(f_left_minus_eps, 
                              f_left_plus_eps, 
                              velocity
                             )
    
    left_flux  = velocity * f_left
    right_flux = np.roll(left_flux, -1)
    
    df_dt = -(right_flux - left_flux) / dr
    return df_dt

t_final = 200
for T in range(1, t_final):
    print('Computing For Time =', T)
    N = 10
    if(T>1):
        while(sol[1]['message'] != 'Integration successful.'):
            print('N is now =', N)
            sol = odeint(df_dt, old_sol, np.linspace(T-1, T, N), full_output = 1)
            N  *= 10
        old_sol = sol[0][-1]
        # Reseting the message
        sol[1]['message'] = '1'
    else:
        sol = odeint(df_dt, f0, np.linspace(T-1, T, N), full_output = 1)
        assert(sol[1]['message'] == 'Integration successful.')
        old_sol = sol[0][-1]
        # Reseting the message
        sol[1]['message'] = '1'
    
    h5f = h5py.File('fvm_data/%04d'%T + '.h5', 'w')
    h5f.create_dataset('time', data = np.linspace(T-1, T, N))
    h5f.create_dataset('sol', data = sol[0])
    h5f.close()
