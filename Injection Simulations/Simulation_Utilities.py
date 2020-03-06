import numpy as NP

def minmod(x, y, z):

    min_of_all = NP.minimum(NP.minimum(abs(x), abs(y)), abs(z))

    signx = NP.sign(x)
    signy = NP.sign(y)
    signz = NP.sign(z)
    
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
    f_i_plus_one  = NP.roll(input_array, -1)
    f_i_minus_one = NP.roll(input_array,  1)
  
    forward_diff  = (f_i_plus_one - input_array  )
    backward_diff = (input_array  - f_i_minus_one)
    central_diff  = backward_diff + forward_diff

    # Sort of user controlled:
    slope_lim_theta = 2

    left   = slope_lim_theta * backward_diff
    center = 0.5 * central_diff
    right  = slope_lim_theta * forward_diff

    return(minmod(left, center, right))

def reconstructMinmod(variable_to_reconstruct):
    
    slope = slope_minmod(variable_to_reconstruct)

    left_face_value  = variable_to_reconstruct - 0.5 * slope
    right_face_value = variable_to_reconstruct + 0.5 * slope
   
    return(left_face_value, right_face_value)

def riemannSolver(left_state, right_state, velocity):
    upwind_state = NP.select([velocity >= 0, velocity < 0], 
                             [left_state, right_state]
                            )
    return upwind_state
