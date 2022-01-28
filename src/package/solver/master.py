import numpy as np
from build_problem import build
from matrices import G_L
###############################################################################
""" 
[] have main take inputs from YAML 
[] make sigma_t and sigma_s dependent on space

"""
###############################################################################

def main():
    
    tfinal = 1.0
    angles = [64]
    Ms = [4]
    N_spaces = [2,4,8,16]
    x0 = 1e-12
    
    source = "plane"
    uncollided = True
    moving = True
    time = True 
    plotting = True
    
    for count, N_space in enumerate(N_spaces):
        sigma_t = np.ones(N_space)
        sigma_s = np.ones(N_space)
        M = Ms[0]
        N_ang = angles[0]
        
        initialize = build(N_ang, N_space, M, tfinal, x0, sigma_t, sigma_s, source, uncollided, moving, time, plotting)
        initialize.make_IC()
        IC = initialize.IC
        matrices = G_L(initialize)
        
main()
        
        
        
    
    
    
    

            





