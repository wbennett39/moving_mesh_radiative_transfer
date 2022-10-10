from numba import njit, types, prange
import quadpy
import ctypes
from numba.extending import get_cython_function_address
import numpy as np
import math
from scipy.special import expi
import matplotlib.pyplot as plt
from ..plots.plot_functions.show import show


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
@njit('float64[:](int64, float64[:], float64, float64)')
def normPn(n,x,a=-1.0,b=1.0):
    tmp = 0*x
    for count in prange(x.size):
        z = (b+a-2*x[count])/(a-b)
        fact = np.sqrt((2*n+1)/(b-a)) #*(x>=a)*(x<=b)
        # tmp[count] = sc.eval_legendre(n,z)*fact
        tmp[count] = numba_eval_legendre_float64(n, z)
    return tmp * fact

@njit
def dx_normPn(n, x, a = -1.0, b = 1.0):
    tmp = 0*x
    fact = np.sqrt((2*n+1)/(b-a))
    for count in prange(x.size):
        z = (b+a-2*x[count])/(a-b)
        if n == 0:
            tmp[count] = 0.0
        elif n == 1:
            tmp[count] = 1.0
        elif n == 2:
            tmp[count] = 3*z
        elif n == 3:
            tmp[count] = (-3 + 15*z**2)/2.
        elif n == 4:
            tmp[count] = (-60*z + 140*z**3)/8.
        elif n == 5:
            tmp[count] = (15 - 210*z**2 + 315*z**4)/8.
        elif n == 6:
            tmp[count] = (210*z - 1260*z**3 + 1386*z**5)/16.
        elif n == 7:
            tmp[count] = (-35 + 945*z**2 - 3465*z**4 + 3003*z**6)/16.
        elif n == 8:
            tmp[count] = (-2520*z + 27720*z**3 - 72072*z**5 + 51480*z**7)/128.
        elif n == 9:
            tmp[count] = (315 - 13860*z**2 + 90090*z**4 - 180180*z**6 + 109395*z**8)/128.
        elif n == 10:
            tmp[count] = (6930*z - 120120*z**3 + 540540*z**5 - 875160*z**7 + 461890*z**9)/256.
        elif n == 11:
            tmp[count] = (-693 + 45045*z**2 - 450450*z**4 + 1531530*z**6 - 2078505*z**8 + 969969*z**10)/256.


    return tmp * fact




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
    if t <= t0:
        if (abx > x0):
            arg1 = max(0,edge - sqrt3 * (abx - x0))
            arg2 = max(0,edge - sqrt3 * (abx + x0))
            # arg2 = max(arg2, edge)
            # arg2 = min(arg2,t0)
            # arg1 = min(arg2,t0)
            return s2_F(t, arg1) - s2_F(t, arg2)
        
        elif (abx <= x0):
            if (edge + sqrt3 * abx <= sqrt3 * x0):
                return  2 * (s2_F(t, edge) - s2_F(t, 0))
            elif (edge + sqrt3 * abx > sqrt3 * x0) and (edge - sqrt3 * (abx + x0) > 0):
                arg2 = edge - sqrt3 * (x0 - abx)
                arg1 = edge - sqrt3 * (abx + x0)
                if arg1 <0 or arg2 <0:
                    print("error negative bounds")
                return s2_F(t, arg2) - s2_F(t, arg1) + 2 * (s2_F(t, edge) - s2_F(t, arg2))
            elif (edge + sqrt3 * abx > sqrt3 * x0) and (edge - sqrt3 * (abx + x0) <= 0): 
                arg1 = max(0,edge - sqrt3 * (x0 - abx))
                if arg1 <0:
                    print("error negative bounds")
                return 2 * (s2_F(t, edge) - s2_F(t, arg1)) + s2_F(t, arg1) - s2_F(t, 0)
            else:
                print("missed case")
    elif t > t0:
        T0 = edge
        x = abs(x)
        if x0 - math.sqrt(3)*(t-T0)/3.0 <= x <= x0 + math.sqrt(3)*(t-T0)/3.0:
            arg = max(T0 - (x-x0)*3/math.sqrt(3),0)
            return s2_F(t,  T0) - s2_F(t, 0) 
            
        elif x > x0 + math.sqrt(3)*(t-T0)/3.0:
            arg = max(t - (x-x0)*3/math.sqrt(3),0)
            return s2_F(t,  arg) - s2_F(t, 0) 
            
        elif x < x0 - math.sqrt(3)*(t-T0)/3.0:
            arg = t - (x0-x)*3/math.sqrt(3)
            if t - (x0-x)*3/math.sqrt(3) <= 0:
                return  2*(s2_F(t,  T0) - s2_F(t, 0))
            elif t - (x0-x)*3/math.sqrt(3) > 0 and t - (x0-x)*3/math.sqrt(3) < T0:
                return 2*(s2_F(t,  T0) - s2_F(t, arg)) + s2_F(t, arg) - s2_F(t, 0)
            else:
                return 0
         
def su_olson_s2_integrand(tau,x,t,x0,t0):
    return  (np.exp(-t + tau)*(-np.heaviside((-3*x - 3*x0 + math.sqrt(3)*(t - tau))/(t - tau),1) - np.heaviside((3*x - 3*x0 + math.sqrt(3)*(t - tau))/(t - tau),1) + np.heaviside((-3*x + 3*x0 + math.sqrt(3)*(t - tau))/(t - tau),1) + np.heaviside((3*x + 3*x0 + math.sqrt(3)*(t - tau))/(t - tau),1)))/2.
    
            
@njit
def uncollided_s2_gaussian(x,t,sigma,t0):
    tf = min(t,t0)

    return (math.exp(-(math.sqrt(3)*x) + (3*sigma**2)/4.)*math.sqrt(3*math.pi)*sigma*(math.erf((-2*t + 2*tf + 2*math.sqrt(3)*x - 3*sigma**2)/(2.*math.sqrt(3)*sigma)) + math.exp(2*math.sqrt(3)*x)*(math.erf(t/(math.sqrt(3)*sigma) + x/sigma + (math.sqrt(3)*sigma)/2.) - math.erf((2*t - 2*tf + 2*math.sqrt(3)*x + 3*sigma**2)/(2.*math.sqrt(3)*sigma))) + math.erf((2*math.sqrt(3)*t - 6*x + 3*math.sqrt(3)*sigma**2)/(6.*sigma))))/4.

    
@njit
def uncollided_s2_gaussian_thick(x,t,sigma,t0):
    return (6*t**5 + t**6 + 12*t**3*(10 + 15*x**2 - 3*sigma**2) + 3*t**4*(10 + 15*x**2 - 3*sigma**2) +  18*t*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) +  9*t**2*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) - 9*(-80 - 3*x**2*(40 + 10*x**2 + x**4) + 24*sigma**2 + 9*x**2*(4 + x**2)*sigma**2 - 6*(2 + 3*x**2)*sigma**4 + 18*sigma**6) + math.exp(t0)*(-t**6 + 6*t**5*(-1 + t0) + 6*t0**5 - t0**6 + 12*t0**3*(10 + 15*x**2 - 3*sigma**2) - 3*t0**4*(10 + 15*x**2 - 3*sigma**2) - 3*t**4*(10 + 5*(-2 + t0)*t0 + 15*x**2 - 3*sigma**2) + 4*t**3*(-30 + 5*t0*(6 + (-3 + t0)*t0) + 45*(-1 + t0)*x**2 - 9*(-1 + t0)*sigma**2) + 18*t0*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) - 9*t0**2*(40 + 15*x**2*(4 + x**2) - 6*(2 + 3*x**2)*sigma**2 + 6*sigma**4) + 9*(-80 - 3*x**2*(40 + 10*x**2 + x**4) + 24*sigma**2 + 9*x**2*(4 + x**2)*sigma**2 - 6*(2 + 3*x**2)*sigma**4 + 18*sigma**6) - 3*t**2*(5*(24 + t0*(-24 + t0*(12 + (-4 + t0)*t0))) + 45*x**4 - 18*(2 + (-2 + t0)*t0)*sigma**2 + 18*sigma**4 + 18*x**2*(10 + 5*(-2 + t0)*t0 - 3*sigma**2)) + 6*t*(-120 + t0*(120 + t0*(-60 + t0*(20 + (-5 + t0)*t0))) + 45*(-1 + t0)*x**4 - 6*(-6 + t0*(6 + (-3 + t0)*t0))*sigma**2 + 18*(-1 + t0)*sigma**4 + 6*x**2*(-30 + 5*t0*(6 + (-3 + t0)*t0) - 9*(-1 + t0)*sigma**2))))/(162.*math.exp(t)*sigma**6)

@njit        
def problem_identifier(source_type):
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

@njit
def uncollided_su_olson_s2_2(x,t,x0,t0):
    T0 = min(t, t0)
    if t > t0:
        x = abs(x)
        if x0 - math.sqrt(3)*(t-T0)/3.0 <= abs(x) <= x0 + math.sqrt(3)*(t-T0)/3.0:
            arg = max(T0 - (abs(x)-x0)*3/math.sqrt(3),0)
            return s2_F(t,  T0) - s2_F(t, 0) 
        
        
        if abs(x) < x0 - math.sqrt(3)*(t-T0)/3.0:
            arg = t - (x0-x)*3/math.sqrt(3)
            
            if t - (x0-x)*3/math.sqrt(3) <= 0:
                return  2*(s2_F(t,  T0) - s2_F(t, 0))
            
            
            elif t - (x0-x)*3/math.sqrt(3) > 0 and t - (x0-x)*3/math.sqrt(3) < T0:
                return 2*(s2_F(t,  T0) - s2_F(t, arg)) + s2_F(t, arg) - s2_F(t, 0)
            
            

            
        elif x > x0 + math.sqrt(3)*(t-T0)/3.0:
            arg = max(t - (x-x0)*3/math.sqrt(3),0)
            return s2_F(t,  arg) - s2_F(t, 0) 
       
        
    else:
        return 0.0

def test_s2_sol(t = 10, t0 = 10):
    import scipy.integrate as integrate
    
    xs = np.linspace(380, 420, 1000)
    phi = xs*0
    phi_test = xs*0
    phi_exact = xs*0
    x0 = 400 
    for ix in range(xs.size):
        phi[ix] = uncollided_su_olson_s2(xs[ix],t, x0, t0)
        phi_test[ix] = uncollided_su_olson_s2_2(xs[ix],t, x0, t0)
        phi_exact[ix] = integrate.quad(su_olson_s2_integrand, 0, min(t,t0), args = (xs[ix],t,x0,t0))[0]
    
    # plt.plot(xs, phi, '-ob')
    plt.plot(xs, phi_exact, '-k')
    # plt.plot(xs, phi_test, '-or', mfc = 'none')
    
    print(np.sqrt(np.mean(phi_exact-phi)**2), 'RMSE')
    show('uncollided_su_olson_s2_t_10')
    plt.show()

def test_square_sol(t = 240, t0 = 10):
    import scipy.integrate as integrate
    
    xs = np.linspace(0, 480, 1000)
    phi = xs*0
    phi_test = xs*0
    phi_exact = xs*0
    x0 = 400 
    for ix in range(xs.size):
        phi[ix] = uncollided_square_s2(xs[ix],t, x0, t0)

    
    # plt.plot(xs, phi, '-ob')
    plt.plot(xs, phi, '-k')
    # plt.plot(xs, phi_test, '-or', mfc = 'none')
    
    show('uncollided_square_s2t')
    plt.show()

# def time_step_counter(t, division):
    
    
    
        
     
        
        
        
        
        
        
        
        
        
        
        
        

