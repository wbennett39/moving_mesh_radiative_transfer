from numba import njit, jit, int64, float64
from numba.experimental import jitclass
import numpy as np
import math


data = [('datax', float64[:]),
        ('datay', float64[:]),
        ('n', int64),
        ('coeff_array', float64[:,:])
        ]
# @jitclass(data) 
class cubic_spline(object):
    def __init__(self, datax, datay):
        self.datax = datax
        self.datay = datay
        self.n = datay.size - 1
        self.coeff_array = np.zeros((4*self.n, 4*self.n))
        self.rhs = np.zeros(4*self.n)
        
        self.solve_coefficients()


    def solve_coefficients(self):  
        rhs = np.zeros(4*self.n)
        #set up matrix equations
        self.points()
        self.first_deriv()
        self.second_deriv()
        
        # first derivative equations 

    def points(self):
        # first row
        x = self.datay[0]
        self.coeff_array[0, 0:4] = [1, x, x**2, x**3]
        self.rhs[0] = self.datay[0]
        count_1 = 1 
        count_2 = 0
        count_3 = 1
        left_index = 0 # indices for the matrix
        right_index = 4
        for ix in range(1, 2*self.n):
            x = self.datax[count_1]
            self.coeff_array[ix, left_index:right_index] = [1, x, x**2, x**3]
            self.rhs[ix] = self.datay[count_1]
            count_2 += 1
            count_3 += 1
            if count_2 == 2:
                count_2 = 0
                count_1 += 1
            if count_3 == 2:
                left_index = right_index
                right_index += 4
    
    def first_deriv(self):
        count_1 = 2*self.n
        count_2 = 1
        left_index = 0
        right_index = 4
        for ix in range(1, self.n):
            x = self.datax[ix]
            self.coeff_array[count_1, left_index:right_index] = [0, 1, 2*x, 3*x**2]
            self.coeff_array[count_1, left_index+4:right_index+4] = [0, -1, -2*x, -3*x**2]
            self.rhs[count_1] = 0
            count_1 += 1
            count_2 += 1
            # if count_2 == 2:
            left_index += 4
            right_index += 4
                # count_2 = 0
    
    def second_deriv(self):
        count_1 = 2*self.n + self.n-2 + 1
        count_2 = 1
        left_index = 0
        right_index = 4
        for ix in range(1, self.n):
            x = self.datax[ix]
            self.coeff_array[count_1, left_index:right_index] = [0, 0, 2, 6*x]
            self.rhs[count_1] = 0
            self.coeff_array[count_1, left_index+4:right_index+4] = [0, 0, -2, -6*x]
            count_1 += 1
            count_2 += 1
            if count_2 == 2:
                left_index += 4
                right_index += 4
                count_2 = 0
        # set first point
        x = self.datax[0]
        self.coeff_array[count_1,0:4 ] = [0, 0, -2, -6*x]
        self.rhs[count_1] = 0
        x = self.datax[-1]
        self.coeff_array[count_1+1,left_index+4:right_index+4] = [0, 0, 2, 6*x]
        self.rhs[count_1+1]=0

        



def test_spline():
    pi = math.pi
    x = np.array([0, pi/2, pi, 3*pi/2, 2*pi ])
    y = np.sin(x)
    spline_object = cubic_spline(x, y)
    rmc_coef, rmc_rhs = mcclarren_spline()
    print(spline_object.coeff_array - rmc_coef)
    np.testing.assert_allclose(rmc_rhs,spline_object.rhs, atol = 1e-15)
    np.testing.assert_allclose(rmc_coef, spline_object.coeff_array, atol = 1e-15 )


    #knot points are sin(x) at 0, pi/2,pi, 3pi/2, 2pi

def mcclarren_spline():
    n = 4 #four intervals
    data = np.array([(0,0),(np.pi*0.5,1),(np.pi,0),(np.pi*1.5,-1),(np.pi*2,0)])
    coef_matrix = np.zeros((4*n,4*n))
    rhs = np.zeros(4*n)
    #set up the 2n equations that match the data at the knot points
    #first point
    x = data[0,0]
    coef_matrix[0,0:4] = [1,x,x**2,x**3]
    rhs[0] = data[0,1]
    #second point
    x = data[1,0]
    coef_matrix[1,0:4] = [1,x,x**2,x**3]
    rhs[1] = data[1,1]
    x = data[1,0]
    coef_matrix[2,4:8] = [1,x,x**2,x**3]
    rhs[2] = data[1,1]
    #third point
    x = data[2,0]
    coef_matrix[3,4:8] = [1,x,x**2,x**3]
    rhs[3] = data[2,1]
    x = data[2,0]
    coef_matrix[4,8:12] = [1,x,x**2,x**3]
    rhs[4] = data[2,1]
    #fourth point
    x = data[3,0]
    coef_matrix[5,8:12] = [1,x,x**2,x**3]
    rhs[5] = data[3,1]
    x = data[3,0]
    coef_matrix[6,12:16] = [1,x,x**2,x**3]
    rhs[6] = data[3,1]
    #last point
    x = data[4,0]
    coef_matrix[7,12:16] = [1,x,x**2,x**3]
    rhs[7] = data[4,1]
    #now the first derivative equations
    #second point
    x = data[1,0]
    coef_matrix[8,0:4] = [0,1,2*x,3*x**2]
    rhs[8] = 0
    coef_matrix[8,4:8] = [0,-1,-2*x,-3*x**2]
    #third point
    x = data[2,0]
    coef_matrix[9,4:8] = [0,1,2*x,3*x**2]
    rhs[9] = 0
    coef_matrix[9,8:12] = [0,-1,-2*x,-3*x**2]
    #fourth point
    x = data[3,0]
    coef_matrix[10,8:12] = [0,1,2*x,3*x**2]
    rhs[10] = 0
    coef_matrix[10,12:16] = [0,-1,-2*x,-3*x**2]
    #now the second derivative equations
    #second point
    x = data[1,0]
    coef_matrix[11,0:4] = [0,0,2,6*x]
    rhs[11] = 0
    coef_matrix[11,4:8] = [0,0,-2,-6*x]
    #third point
    x = data[2,0]
    coef_matrix[12,4:8] = [0,0,2,6*x]
    rhs[12] = 0
    coef_matrix[12,8:12] = [0,0,-2,-6*x]
    #fourth point
    x = data[3,0]
    coef_matrix[13,4:8] = [0,0,2,6*x]
    rhs[13] = 0
    coef_matrix[13,8:12] = [0,0,-2,-6*x]
    #set first point to 0
    x = data[0,0]
    coef_matrix[14,0:4] = [0,0,-2,6*x]
    rhs[14] = 0
    #set last point to 0
    x = data[4,0]
    coef_matrix[15,12:16] = [0,0,2,6*x]
    rhs[15] = 0
    return coef_matrix, rhs


                


test_spline()