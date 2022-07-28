from numba import njit, types, prange
import quadpy
import ctypes
from numba.extending import get_cython_function_address
import numpy as np
import math
from scipy.special import expi
import matplotlib.pyplot as plt


@njit 
def integrate_quad(a, b, xs, ws, func1, func2):
    return (b-a)/2 * np.sum(ws * func1((b-a)/2*xs + (a+b)/2) * func2((b-a)/2*xs + (a+b)/2))

_dble = ctypes.c_double
addr = get_cython_function_address("scipy.special.cython_special", "__pyx_fuse_0_1eval_legendre")
functype = ctypes.CFUNCTYPE(_dble, _dble, _dble)
eval_legendre_float64_fn = functype(addr)

# @njit("float64[:](float64,float64[:])")  
@njit
def numba_eval_legendre_float64(n, x):
      return eval_legendre_float64_fn(n, x)
  
addr = get_cython_function_address("scipy.special.cython_special", "__pyx_fuse_1expi")
functype = ctypes.CFUNCTYPE(_dble, _dble)
expn_fn = functype(addr)

@njit("float64(float64)")
def numba_expi(x):
    return expn_fn(x)

# @njit("float64[:](float64,float64[:],float64,float64)", looplift=False, parallel=False)
@njit
def normPn(n,x,a=-1.0,b=1.0):
    tmp = 0*x
    for count in prange(x.size):
        z = (b+a-2*x[count])/(a-b)
        fact = np.sqrt((2*n+1)/(b-a)) #*(x>=a)*(x<=b)
        # tmp[count] = sc.eval_legendre(n,z)*fact
        tmp[count] = numba_eval_legendre_float64(n, z)*fact
    return tmp
    

# @njit("float64[:](int64, float64[:], float64[:], float64[:,:,:], int64, float64[:,:])", parallel = True, looplift = True, fastmath = True)
def make_phi(N_ang, ws, xs, u, M, edges):
    output = xs*0
    psi = np.zeros((N_ang, xs.size))
    for ang in range(N_ang):
        for count in range(xs.size):
            idx = np.searchsorted(edges[:], xs[count])
            if (idx == 0):
                idx = 1
            if (idx >= edges.size):
                idx = edges.size - 1
            for i in range(M+1):
                psi[ang, count] += u[ang,idx-1,i] * normPn(i,xs[count:count+1],float(edges[idx-1]),float(edges[idx]))[0]
    output = np.sum(np.multiply(psi.transpose(), ws), axis = 1)
    return output
@njit 
def surf_func(speed,u,space,j,side,xL,xR,N_space):
#    print(side)
    if j ==0:
        B_right = 1/math.sqrt(xR-xL)
        B_left = 1/math.sqrt(xR-xL)
    else:
         B_right = math.sqrt(2*j+1)/math.sqrt(xR-xL)
         if j%2 ==0:
                B_left = math.sqrt(2*j+1)/math.sqrt(xR-xL)
         else:
                B_left = -math.sqrt(2*j+1)/math.sqrt(xR-xL)
    if speed == 0:
        return 0
    elif speed > 0 and side == "R":
        return u[space,j]*B_right
    elif speed > 0 and side =="L":
        if space !=0:
#            print(u[k-1,j])
            return u[space-1,j]*B_right 
        else:
            return 0
    elif speed < 0 and side =="R":
        if space != N_space-1:
            return u[space+1,j]*B_left
        else:
            return 0
    elif speed < 0 and side =="L":
        return u[space,j]*B_left

@njit
def LU_surf_func(u,space,N_space,mul,M,xL,xR,dxL,dxR):
    sumright = 0
    sumleft = 0
    rightspeed = mul - dxR
    leftspeed = mul-dxL
    for j in range(0,M+1):
        sumright += surf_func(rightspeed,u,space,j,"R",xL,xR,N_space)
        sumleft += surf_func(leftspeed,u,space,j,"L",xL,xR,N_space)
    LU = np.zeros(M+1).transpose()
    for i in range(0,M+1):
        if i == 0:
            B_right = 1/math.sqrt(xR-xL)
            B_left = 1/math.sqrt(xR-xL)
        elif j>0:
            B_right = math.sqrt(2*i+1)/math.sqrt(xR-xL)
            if i%2 ==0:
                B_left = math.sqrt(2*i+1)/math.sqrt(xR-xL)
            else: 
                B_left = -math.sqrt(2*i+1)/math.sqrt(xR-xL)
        LU[i] = rightspeed*B_right*(sumright) - leftspeed*B_left*(sumleft)
    return LU 

def find_nodes(edges, M):
    scheme = quadpy.c1.gauss_legendre(M+1)
    xs_quad = scheme.points
    ixx = xs_quad.size
    xs_list = np.zeros(ixx*(edges.size-1))
    for k in range(edges.size-1):
        xL = edges[k]
        xR = edges[k+1]
        xs_list[k*ixx:(k+1)*ixx] = xL + (xs_quad + 1)*(xR - xL)/2
    return xs_list

def convergence(err1, x1, err2, x2):
    return -math.log(err2/err1)/math.log(x2/x1)
@njit
def f1(t, tau, x0):
    return -x0 * numba_expi(tau-t)
@njit    
def f2(t, tau, x0, x):
    if tau != t:
        return 0.5*((-x0 + abs(x)) * numba_expi(tau-t) + math.exp(tau - t))
    else:
        return 0.5 
@njit
def f3(t, tau, x0, x):
    return math.exp(tau-t)
@njit
def uncollided_square_s2(x, t, x0, t0):
    t_ceiling = min(t,t0)
    if t > 0:
        tau_1 = 0.0
        end = min(t_ceiling, t - abs(x) + x0)
        if end <= 0.0:
            return 0.0
        tau_2 = min(end, t - x0 - abs(x))
        if tau_2 < 0.0:
            tau_2 = 0.0
        tau_3 = min(end, t - x0 + abs(x))
        if tau_3 < 0.0:
            tau_3 = 0.0
        tau_4 = end
        if tau_4 < 0.0:
            tau_4 = 0.0
        t1 = f1(t, tau_2, x0) - f1(t, tau_1, x0)
        t2 = f2(t, tau_3, x0, x) - f2(t, tau_2, x0, x)
        t3 = f3(t, tau_4, x0, x) - f3(t, tau_3, x0, x)
        
        return t1 + t2 + t3
    else:
        return 0.0

@njit
def s2_F(t,tau):
    """ integrand for uncollided square s2
    """
    return 0.5 * math.exp(-t + tau)

@njit 
def uncollided_su_olson_s2(x,t,x0,t0):
    sqrt3 = math.sqrt(3)
    abx = abs(x)
    edge = min(t,t0)
    if (abx > x0):
        arg1 = max(0,t - sqrt3 * (abx - x0))
        arg2 = max(0,t - sqrt3 * (abx + x0))
        return s2_F(t, arg1) - s2_F(t, arg2)
    
    elif (abx <= x0):
        if (t + sqrt3 * abx <= sqrt3 * x0):
            return  2 * (s2_F(t, edge) - s2_F(t, 0))
        elif (t + sqrt3 * abx > sqrt3 * x0) and (t - sqrt3 * (abx + x0) > 0):
            arg2 = t - sqrt3 * (x0 - abx)
            arg1 = t - sqrt3 * (abx + x0)
            if arg1 <0 or arg2 <0:
                print("error negative bounds")
            return s2_F(t, arg2) - s2_F(t, arg1) + 2 * (s2_F(t, edge) - s2_F(t, arg2))
        elif (t + sqrt3 * abx > sqrt3 * x0) and (t - sqrt3 * (abx + x0) <= 0): 
            arg1 = max(0,t - sqrt3 * (x0 - abx))
            if arg1 <0:
                print("error negative bounds")
            return 2 * (s2_F(t, edge) - s2_F(t, arg1)) + s2_F(t, arg1) - s2_F(t, 0)
        else:
            print("missed case")
            
@njit
def uncollided_s2_gaussian(x,t,sigma,t0):
    tf = min(t,t0)

    return (math.exp(-(math.sqrt(3)*x) + (3*sigma**2)/4.)*math.sqrt(3*math.pi)*sigma*(math.erf((-2*t + 2*tf + 2*math.sqrt(3)*x - 3*sigma**2)/(2.*math.sqrt(3)*sigma)) + math.exp(2*math.sqrt(3)*x)*(math.erf(t/(math.sqrt(3)*sigma) + x/sigma + (math.sqrt(3)*sigma)/2.) - math.erf((2*t - 2*tf + 2*math.sqrt(3)*x + 3*sigma**2)/(2.*math.sqrt(3)*sigma))) + math.erf((2*math.sqrt(3)*t - 6*x + 3*math.sqrt(3)*sigma**2)/(6.*sigma))))/4.

    
@njit
def uncollided_s2_gaussian_thick(x,t,sigma,t0):
    return (6*t**5 + t**6 + 12*t**3*(10 + 15*x**2 - 3*sigma**2) + 3*t**4*(10 + 15*x**2 - 3*sigma**2) +  18*t*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) +  9*t**2*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) - 9*(-80 - 3*x**2*(40 + 10*x**2 + x**4) + 24*sigma**2 + 9*x**2*(4 + x**2)*sigma**2 - 6*(2 + 3*x**2)*sigma**4 + 18*sigma**6) + math.exp(t0)*(-t**6 + 6*t**5*(-1 + t0) + 6*t0**5 - t0**6 + 12*t0**3*(10 + 15*x**2 - 3*sigma**2) - 3*t0**4*(10 + 15*x**2 - 3*sigma**2) - 3*t**4*(10 + 5*(-2 + t0)*t0 + 15*x**2 - 3*sigma**2) + 4*t**3*(-30 + 5*t0*(6 + (-3 + t0)*t0) + 45*(-1 + t0)*x**2 - 9*(-1 + t0)*sigma**2) + 18*t0*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) - 9*t0**2*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) + 9*(-80 - 3*x**2*(40 + 10*x**2 + x**4) + 24*sigma**2 + 9*x**2*(4 + x**2)*sigma**2 - 6*(2 + 3*x**2)*sigma**4 + 18*sigma**6) - 3*t**2*(5*(24 + t0*(-24 + t0*(12 + (-4 + t0)*t0))) + 45*x**4 - 18*(2 + (-2 + t0)*t0)*sigma**2 + 18*sigma**4 + 18*x**2*(10 + 5*(-2 + t0)*t0 - 3*sigma**2)) + 6*t*(-120 + t0*(120 + t0*(-60 + t0*(20 + (-5 + t0)*t0))) + 45*(-1 + t0)*x**4 - 6*(-6 + t0*(6 + (-3 + t0)*t0))*sigma**2 + 18*(-1 + t0)*sigma**4 + 6*x**2*(-30 + 5*t0*(6 + (-3 + t0)*t0) - 9*(-1 + t0)*sigma**2))))/(162.*math.exp(t)*sigma**6)

@njit        
def problem_identifier(source_array):
        if source_type[0] == 1:
        problem_type = 'plane_IC'
    elif source_type[1] == 1:
        problem_type = 'square_IC'
    elif source_type[2] == 1:
        problem_type = 'square_source'
    elif source_type[3] == 1:
        problem_type = 'gaussian_IC'
    elif source_type[4] == 1:
        problem_type = 'gaussian_source'
    else:
        problem_type =='none'
    return problem_type


# xs  = np.linspace(-1.5,1.5, 100)
# phi = xs*0
# phi_u = xs*0
# tf = 1.0
# expi1 = xs*0
# expi2 = xs*0

# for i in range(len(xs)):
#     expi1[i] = numba_expi(xs[i])
#     expi2[i] = expi(xs[i])

#     # phi_u[i] = uncollided_square_s2(xs[i], tf, 0.5, tf)
# plt.figure(1)

# plt.plot(xs, expi1, "--")
# plt.plot(xs, expi2, ":")

