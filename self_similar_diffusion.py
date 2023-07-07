from moving_mesh_transport.loading_and_saving.load_solution import load_sol
import matplotlib.pyplot as plt
import numpy as np
from marshak_point_source import exact_sol
import math


loader_transport = load_sol('transport', 'plane_IC', 'transport', s2 = True, file_name = 'run_data.hdf5')
loader = load_sol('su_olson_thick_s2', 'plane_IC', 'transfer', s2 = True, file_name = 'run_data.hdf5')

tfinal_list = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0]
RMS_list = []
RMS_list_transport = []
RMS_list_transport_2 = []
tfinal_list_2 = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0]
# tfinal_list_2 = [5000.0]
xi, f, m, k, n, p = exact_sol()

sigma = 800

beta = 4

rho0 = 1
c=29.998
alpha = 0

lambda_ = 0 

mu = 1

f1 = 4 *  0.0137225 # GJ/cm$^3$/keV$^4

a = 0.0137225
annoying_term = lambda_ + 1 + (1-mu)*(alpha+4)/beta 
# A = 16 * sigma / (3 * beta  * f1**((4+alpha)/beta)) * rho0 ** annoying_term
A = c/3/sigma
A2 = c/3/sigma/2
A_transport = 3/c

Q = 1 * c * a * sigma

twokm = 2- k- m

def RMSE(l1, l2):
    diff = (l1-l2)**2
    return np.sqrt(np.mean(diff))

def benchmark(xi):
    return 0.5 * (1/math.sqrt(math.pi)) * np.exp(-xi**2/4)

def ss_suolson(N_spaces):
    for it, t in enumerate(tfinal_list):
        loader.call_sol(t, 4, 1e-14, N_spaces, 'rad', False, False)
        loader_transport.call_sol(t, 4, 1e-14, 64, 'rad', False, False)
        x = loader.xs
        # print(x)
        tau = t /c 
        tau_transport = t / c / sigma
        x_transport = loader_transport.xs / sigma
        xi_sol = x / math.sqrt(A*tau)

        dimensionalize = (1 / (A * tau))**(1/2)
            
        dimensionalize_transport = ((1) / (A * tau_transport))**(1/2)

        scaled_phi = loader.phi 

        transport_phi = loader_transport.phi

        dimensional_energy = scaled_phi 
        scaled_e = loader.e
        # scaled_e = scaled_phi * 4 * a
        nondim = scaled_e / dimensionalize
        RMS_list.append(RMSE(dimensional_energy / dimensionalize, benchmark(xi_sol)))

        xi_sol_transport = x_transport / ((Q**n * (A*tau_transport)) ** (1/p))
        plt.figure(1)
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        if t == tfinal_list[-1]:
            plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none')
        plt.plot(xi_sol, dimensional_energy / dimensionalize   , label = f't = {t}')
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$f(\xi)$')
        plt.xlim(-10,10)
        plt.ylim(0, 0.6)

        plt.figure(2)
        plt.plot(x, dimensional_energy, label = f't = {t}')
        plt.xlabel('x')
        plt.ylabel('e')
    # plt.show()

    plt.legend()
    plt.show()
    plt.figure(5)
    plt.loglog(tfinal_list, RMS_list, '-o')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.show()



def ss_transport(N_spaces):
    RMS_list_transport = []
    for it, t in enumerate(tfinal_list_2):
        loader_transport.call_sol(t, 8, 1e-14, N_spaces, 'rad', False, False)
        x = loader_transport.xs / sigma
        x2 = loader_transport.xs / sigma/2
        # print(x)
        tau = t /c /sigma
        tau2 = t /c /sigma/2
        dimensionalize = (1 / (A * tau))**(1/2)
        dimensionalize2 = (1 / (A2 * tau2))**(1/2)

        transport_phi = loader_transport.phi 
        xi_sol = x / math.sqrt(A*tau)
        xi_sol2 = x2 / math.sqrt(A2*tau2)

        plt.figure(3)
        plt.ion()
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        if t == tfinal_list_2[-1]:
            plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        plt.plot(xi_sol, transport_phi / dimensionalize * sigma   , label = f't = {t}, {N_spaces} spatial cells')
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize * sigma , benchmark(xi_sol)))
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        plt.xlim(-10,10)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$f(\xi)$')


    plt.legend()
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.show()


def ss_special_IC(N_spaces):
    RMS_list_transport = []
    for it, t in enumerate(tfinal_list_2):
        loader_transport.call_sol(t, 8, 1e-1, N_spaces, 'rad', False, False)
        x = loader_transport.xs / sigma
        x2 = loader_transport.xs / sigma/2
        # print(x)
        tau = t /c /sigma
        tau2 = t /c /sigma/2
        dimensionalize = (1 / (A * tau))**(1/2)
        dimensionalize2 = (1 / (A2 * tau2))**(1/2)

        transport_phi = loader_transport.phi 
        xi_sol = x / math.sqrt(A*tau)
        xi_sol2 = x2 / math.sqrt(A2*tau2)

        plt.figure(3)
        plt.ion()
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        if t == tfinal_list_2[-1]:
            plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        plt.plot(xi_sol, transport_phi / dimensionalize * sigma   , label = f't = {t}, {N_spaces} spatial cells')
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize * sigma , benchmark(xi_sol)))
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        plt.xlim(-10,10)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$f(\xi)$')


    plt.legend()
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.show()

def ss_dipole(N_spaces):
    RMS_list_transport = []
    for it, t in enumerate(tfinal_list_2):
        loader_transport.call_sol(t, 8, 0.25, N_spaces, 'rad', False, False)
        x = loader_transport.xs / sigma
        x2 = loader_transport.xs / sigma/2
        # print(x)
        tau = t /c /sigma
        tau2 = t /c /sigma/2
        x0 = 0.25
        t0 = -x0**2/6/A
        # t0 = 0
        # dimensionalize = x0*np.abs(x-x0) / ((A * (tau-t0))**1.5)
        dimensionalize = 1/(A * (tau-t0) /x0/800)


        print(x0/A/(tau-t0)*0.25, 'max phi predicted')
        
        print(x0**2/A, 'reasonable time')
        transport_phi = loader_transport.phi 
        xi_sol = (x) / math.sqrt(A*(tau-t0))
        print(max(transport_phi), 'max phi code')

        plt.figure(3)
        plt.ion()
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        # if t == tfinal_list_2[-1]:
        #     plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        plt.plot(xi_sol, transport_phi / dimensionalize * 800     , label = f't = {t}, {N_spaces} spatial cells')
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize  * 800, benchmark(xi_sol)))
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        plt.xlim(0,10)
        plt.ylim(0,1)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$f(\xi)$')


    plt.legend()
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.show()


# ss_transport(32)
# ss_transport(64)
# ss_transport(96)
# ss_transport(128)
# ss_transport(192)
# ss_transport(256)
# ss_transport(384)
# ss_transport(512)

# ss_suolson(64)
# ss_suolson(128)

# ss_dipole(64)
ss_dipole(128)







