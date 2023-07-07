from numba import njit, types, prange
from numba.extending import get_cython_function_address
import numpy as np
import math
from numba import guvectorize
import scipy.interpolate



@njit
def positivize(wavespeed_array):
    # for it in range(1, wavespeed_array.size):
    #     if wavespeed_array[it] 
    return 0.0
        


def set_func(self, x_index, loc, speed):
    """
    This function takes an array of indexes and initializes the 
    edges and Dedges at the indexes given to the respective speeds and locations
    """
    for count, index in enumerate(x_index):
        self.edges[int(index)] = loc[count]
        self.Dedges[int(index)] = speed[count]

def thick_square_moving_func(self, t):
    """ Depreciated function
    """
    # print(self.edges)
    # print(t, "t")
    # print("-------------------------------------")
    half = int(self.N_space/2)
    index = np.searchsorted(self.wave_loc_array[0,0,:], t)
    # self.delta_t = self.wave_loc_array[0,0,index+1] - t 
    if self.tactual == 0.0:
        delim = '#########################################################'
        # print(delim)
        # print(self.wave_loc_array[0,2,:])
        # print(self.wave_loc_array[0,1,:])
        # print(delim)

        if t == 1.0:
            pad = 0.5
        elif t == 100.0:
            pad = 2.0
        else:
            pad = 0.0
        self.delta_t = self.tfinal
        index -=1
        self.right_speed = (-self.wave_loc_array[0,2,index+1] - pad - self.edges[0])/self.delta_t 
        self.left_speed = (-self.wave_loc_array[0,1,index+1]  - self.edges[half-1])/self.delta_t  
        wave_front_array =  self.Dedges_const[0:int(half/2)]*self.right_speed/self.speed*(-1)
        # print(wave_front_array, 'wave speeds forward')
        wave_back_array =  self.Dedges_const[int(half/2)+1:half]*self.left_speed/self.speed
        self.Dedges[0:int(half/2)] = wave_front_array
        self.Dedges[int(half/2)+1:half] = wave_back_array
        self.Dedges[half+1:half+int(half/2)] = - np.flip(wave_back_array)
        self.Dedges[half+int(half/2)+1:] = - np.flip(wave_front_array)

        self.edges = self.edges0 + self.Dedges * self.delta_t
        # print(self.left_speed, 'ls')
        # print(self.wave_loc_array[0,1,index+1], 'left edge')

    else:
        if index != self.index_old or index == 0:
            self.delta_t = self.wave_loc_array[0,0,index+1] - self.wave_loc_array[0,0,index]
            self.right_speed = (self.wave_loc_array[0,2,index+1] - self.wave_loc_array[0,2,index])/self.delta_t
            self.left_speed = (self.wave_loc_array[0,1,index+1] - self.wave_loc_array[0,1,index])/self.delta_t
            # print(self.wave_loc_array[0,2,index+1])
            self.right_speed = self.right_speed*(-1)
            self.left_speed = self.left_speed*(-1)
            
        
        wave_front_array =  self.Dedges_const[0:int(half/2)]*self.right_speed/self.speed*(-1)
            # print(wave_front_array, 'wave speeds forward')
        wave_back_array =  self.Dedges_const[int(half/2)+1:half]*self.left_speed/self.speed
            # print(wave_back_array, 'wave speeds back')
        self.Dedges[0:int(half/2)] = wave_front_array
        self.Dedges[int(half/2)+1:half] = wave_back_array
        self.Dedges[half+1:half+int(half/2)] = - np.flip(wave_back_array)
        self.Dedges[half+int(half/2)+1:] = - np.flip(wave_front_array)

        self.edges = self.edges0 + self.Dedges * (t-self.told)
        self.edges0 = self.edges
        self.told = t
    # print(self.Dedges_const, "const")
    # # old func


@njit
def _interp1d(xnew, xvals, yvals, ynew):
    i = 0
    N = len(xvals)
    if xnew[0] < xvals[0]:
        x_a = 0.0
        y_a = 0.0
        x_b = xvals[0]
        y_b = yvals[0]
    else:
        while xnew[0] >= xvals[i] and i < N:
            i += 1
        if xnew[0] == xvals[i]:
            ynew[0] = yvals[i]
            return ynew
        if i == N:
            i = N-1
        x_a = xvals[i-1]
        y_a = yvals[i-1]
        x_b = xvals[i]
        y_b = yvals[i]
    slope = (xnew[0] - x_a)/(x_b - x_a)
    ynew[0] = slope * (y_b-y_a) + y_a
    return ynew

  
# interp1d_numba = guvectorize(
#     ['float64[:], float64[:], float64[:], float64[:]'],
#     "(),(n),(n) -> ()", nopython=True)(_interp1d)

@njit
def boundary_source_init_func_outside(v0, N_space, x0, tfinal):
        mid = int(N_space/2)
        edges = np.linspace(-x0, x0, N_space+1)
        Dedges = np.copy(edges)*0
        print(v0, 'v0 in init outside')
        # self.Dedges[mid] = - self.fake_sedov_v0
        final_shock_point = - tfinal * v0
        final_edges_left_of_shock = np.linspace(-x0, final_shock_point, int(N_space/2+1))
        final_edges_right_of_shock = np.linspace(final_shock_point, x0, int(N_space/2+1))
        final_edges = np.concatenate((final_edges_left_of_shock[:-1], final_edges_right_of_shock))
        Dedges = (final_edges - edges) / tfinal

        # self.Dedges = self.edges/self.edges[-1] * self.speed 
        
        edges0 = edges

        return edges, edges0, Dedges

@njit
def mesh_cluster_blast_wave(x0, N_spaces, v0, t):
    mid = int(N_space/2)
    shock_loc = v0 * t
    xs = (np.linspace(0, x0, mid + 1))**2
    

