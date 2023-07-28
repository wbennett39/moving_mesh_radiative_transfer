from numba import njit, jit, int64, float64
from numba.experimental import jitclass
import numpy as np
import math
import matplotlib.pyplot as plt


data = [('datax', float64[:]),
        ('datay', float64[:]),
        ('n', int64),
        ('coeff_array', float64[:,:]),
        ('rhs', float64[:]),
        ('coef', float64[:])
        ]
@jitclass(data) 
class cubic_spline(object):
    def __init__(self, datax, datay):
        self.datax = datax
        self.datay = datay
        self.n = datay.size - 1
        print(self.n, 'n')
        self.coeff_array = np.zeros((4*self.n, 4*self.n))
        self.rhs = np.zeros(4*self.n)
        
        self.solve_coefficients()

    def solve_coefficients(self):  
        rhs = np.zeros(4*self.n)
        #set up matrix equations
        self.points()
        self.first_deriv()
        self.second_deriv()
        # compute coefficients
        self.solve_spline()
        # print(self.coeff_array[:, 16:])


    def points(self):
        # first row
        x = self.datax[0]
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
            # print(ix, count_1)
            count_2 += 1
            count_3 += 1
            if count_2 == 2:
                count_2 = 0
                count_1 += 1
            if count_3 == 2:
                left_index += 4
                right_index += 4
                count_3 = 0

    
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
            # count_2 += 1
            # if count_2 == 2:
            left_index += 4
            right_index += 4
            # print(left_index, right_index)
                # count_2 = 0
    
    def second_deriv(self):
        count_1 = 3*self.n-1
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
            
            left_index += 4
            right_index+=4
            if count_2 == 2:
               
                count_2 = 0
        # set first point
        x = self.datax[0]
        self.coeff_array[4*self.n-2, 0:4] = [0, 0, 2, 6*x]
        self.rhs[-2] = 0
        #set last point
        x = self.datax[-1]
        self.coeff_array[4*self.n-1,self.n*4-4:self.n*4] = [0, 0, 2, 6*x]
        self.rhs[-1] = 0
  


    def solve_spline(self):
        # self.coef = np.linalg.solve(self.coeff_array, self.rhs)
        self.coef = GaussElim(self.coeff_array, self.rhs)
        # for it1 in range(self.n*4):
        #     if abs(self.coef[it1] > 1e10):
        #         print(it1, 'it1')
    
    def eval_spline(self, x):
        y_interp = np.zeros(x.size)
        for i in range(x.size):
            for knot in range(self.n):
                if self.datax[knot]<=x[i]<=self.datax[knot+1]:
                    y_interp[i] = np.sum(self.coef[4*knot:(4*knot+4)] * np.array([1,x[i],x[i] **2,x[i]**3]))
                    break #exit for loop
        return y_interp
        

def test_spline():
    pi = math.pi
    x = np.array([0, pi/2, pi, 3*pi/2, 2*pi, 5*pi/2 ])
    x = np.linspace(0, 2*pi, 50)
    y = np.sin(x)
    spline_object = cubic_spline(x, y)
    rmc_coef, rmc_rhs = mcclarren_spline()
    coeffs = GaussElim(rmc_coef, rmc_rhs)
    
    # np.testing.assert_allclose(rmc_rhs,spline_object.rhs, atol = 1e-15)
    # np.testing.assert_allclose(rmc_coef, spline_object.coeff_array, atol = 1e-15 )

    xtest = np.linspace(0, 2*pi, 100)
    ytest = spline_object.eval_spline(xtest)
    err_me = np.sqrt(np.mean((ytest-np.sin(xtest))**2))
    # err_mc = np.abs()
    print(err_me, 'error')

    plt.figure(1)
    plt.plot(xtest, ytest, 'o', mfc = 'none', label = 'interpolated')
    plt.scatter(x, y)
    plt.plot(xtest, np.sin(xtest), label = 'sin(x)')
    plt.legend()
    plt.show()

    plt.figure(2)
    x = np.linspace(-1, 1,500)
    # x = np.array([0.0, 1.0, 2.0])
    # y = np.array([1.0, 3.0, 2.0])
    f1 = lambda x: x**2 * np.exp(-x/4) * np.heaviside(0.5-np.abs(x),0)
    y = f1(x)
    spline_object = cubic_spline(x, y)
    xtest = np.linspace(-1, 1, 100)
    ytest = spline_object.eval_spline(xtest)
    plt.plot(xtest,f1(xtest))
    plt.plot(xtest,ytest, label = 'interpolated')
    plt.scatter(x, y)
    plt.legend()
    plt.show()

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
    coef_matrix[13,8:12] = [0,0,2,6*x]
    rhs[13] = 0
    coef_matrix[13,12:16] = [0,0,-2,-6*x]
    #set first point to 0
    x = data[0,0]
    coef_matrix[14,0:4] = [0,0,2,6*x]
    rhs[14] = 0
    #set last point to 0
    x = data[4,0]
    coef_matrix[15,12:16] = [0,0,2,6*x]
    rhs[15] = 0
    return coef_matrix, rhs

def mcclarren_spline_2():
    n = 5 #four intervals
    data = np.array([(0,0),(np.pi*0.5,1),(np.pi,0),(np.pi*1.5,-1),(np.pi*2,0), (5*np.pi/2, 1.0)])
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
    #fifth point
    x = data[4,0]
    coef_matrix[7,12:16] = [1,x,x**2,x**3]
    rhs[7] = data[4,1]
    x = data[4,0]
    coef_matrix[8,16:20] = [1,x,x**2,x**3]
    rhs[8] = data[4,1]
    #last point
    x = data[5,0]
    coef_matrix[9,16:20] = [1,x,x**2,x**3]
    rhs[9] = data[5,1]

    #now the first derivative equations

    #second point
    x = data[1,0]
    coef_matrix[10,0:4] = [0,1,2*x,3*x**2]
    rhs[10] = 0
    coef_matrix[10,4:8] = [0,-1,-2*x,-3*x**2]
    #third point
    x = data[2,0]
    coef_matrix[11,4:8] = [0,1,2*x,3*x**2]
    rhs[11] = 0
    coef_matrix[11,8:12] = [0,-1,-2*x,-3*x**2]
    #fourth point
    x = data[3,0]
    coef_matrix[12,8:12] = [0,1,2*x,3*x**2]
    rhs[12] = 0
    coef_matrix[12,12:16] = [0,-1,-2*x,-3*x**2]
    #fifth point
    x = data[4,0]
    coef_matrix[13,12:16] = [0,1,2*x,3*x**2]
    rhs[13] = 0
    coef_matrix[13,16:20] = [0,-1,-2*x,-3*x**2]

    #now the second derivative equations

    #second point
    x = data[1,0]
    coef_matrix[14,0:4] = [0,0,2,6*x]
    rhs[14] = 0
    coef_matrix[14,4:8] = [0,0,-2,-6*x]
    #third point
    x = data[2,0]
    coef_matrix[15,4:8] = [0,0,2,6*x]
    rhs[15] = 0
    coef_matrix[15,8:12] = [0,0,-2,-6*x]
    #fourth point
    x = data[3,0]
    coef_matrix[16,4:8] = [0,0,2,6*x]
    rhs[16] = 0
    coef_matrix[16,8:12] = [0,0,-2,-6*x]
    #fifth point
    x = data[4,0]
    coef_matrix[17,8:12] = [0,0,2,6*x]
    rhs[17] = 0
    coef_matrix[17,12:16] = [0,0,-2,-6*x]
    #set first point to 0
    x = data[0,0]
    coef_matrix[18,0:4] = [0,0,2,6*x]
    rhs[18] = 0
    #set last point to 0
    x = data[5,0]
    coef_matrix[19,16:20] = [0,0,2,6*x]
    rhs[19] = 0
    return coef_matrix, rhs


                
@njit
def swap_rows(A, a, b):
    """Rows two rows in a matrix, switch row a with row b
    
    args: 
    A: matrix to perform row swaps on
    a: row index of matrix
    b: row index of matrix
    
    returns: nothing
    
    side effects:
    changes A to rows a and b swapped
    """
    assert (a>=0) and (b>=0)
    N = A.shape[0] #number of rows
    assert (a<N) and (b<N) #less than because 0-based indexing
    temp = A[a,:].copy()
    A[a,:] = A[b,:].copy()
    A[b,:] = temp.copy()
@njit
def BackSub(aug_matrix,x):
    """back substitute a N by N system after Gauss elimination

    Args:
        aug_matrix: augmented matrix with zeros below the diagonal
        x: length N vector to hold solution
    Returns:
        nothing
    Side Effect:
        x now contains solution
    """
    N = x.size
    for row in range(N-1,-1,-1):
        RHS = aug_matrix[row,N]
        for column in range(row+1,N):
            RHS -= x[column]*aug_matrix[row,column]
        x[row] = RHS/aug_matrix[row,row]
    return
@njit
def GaussElim(A,b,LOUD=0):
    """create a Gaussian elimination with pivoting matrix for a system

    Args:
        A: N by N array
        b: array of length N
    Returns:
        solution vector in the original order
    """
    [Nrow, Ncol] = A.shape
    assert Nrow == Ncol
    N = Nrow
    #create augmented matrix
    aug_matrix = np.zeros((N,N+1))
    aug_matrix[0:N,0:N] = A
    aug_matrix[:,N] = b
    #augmented matrix is created
    
    #create scale factors 
    s = np.zeros(N)
    count = 0
    for row in aug_matrix[:,0:N]: #don't include b
        s[count] = np.max(np.fabs(row))
        count += 1
    if LOUD:
        print("s =",s)
    if LOUD:
        print("Original Augmented Matrix is\n",aug_matrix)
    #perform elimination
    for column in range(0,N):
        
        #swap rows if needed
        largest_pos = np.argmax(np.fabs(aug_matrix[column:N,column]/s[column])) + column
        if (largest_pos != column):
            if (LOUD):
                print("Swapping row",column,"with row",largest_pos)
                print("Pre swap\n",aug_matrix)
            swap_rows(aug_matrix,column,largest_pos)
            #re-order s
            tmp = s[column]
            s[column] = s[largest_pos]
            s[largest_pos] = tmp
            if (LOUD):
                print("A =\n",aug_matrix)
        #finish off the row
        for row in range(column+1,N):
            mod_row = aug_matrix[row,:]
            mod_row = mod_row - mod_row[column]/aug_matrix[column,column]*aug_matrix[column,:]
            aug_matrix[row] = mod_row
    #now back solve
    x = b.copy()
    if LOUD:
        print("Final aug_matrix is\n",aug_matrix)
    BackSub(aug_matrix,x)
    return x


test_spline()


