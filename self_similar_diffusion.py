from moving_mesh_transport.loading_and_saving.load_solution import load_sol
import matplotlib.pyplot as plt
import numpy as np
from marshak_point_source import exact_sol
import math
from moving_mesh_transport.plots.plot_functions.show import show

loader_transport = load_sol('transport', 'plane_IC', 'transport', s2 = True, file_name = 'run_data.hdf5')
loader = load_sol('su_olson_thick_s2', 'plane_IC', 'transfer', s2 = True, file_name = 'run_data.hdf5')



def loglog_converge(x, error, start):
    # start = 5
    array_of_slopes = np.zeros(x.size-start-1)
    for ix in range(start, x.size-1):
        # if error[ix-1] > error[ix]:
        t1 = x[ix]
        t2 = x[ix-1]
        y1 = error[ix]
        y2 = error[ix-1]
        slope = (y2-y1)/(t2-t1)
        array_of_slopes[ix - start] = -math.log(y2/y1)/math.log(t2/t1)
        print('######')
        print(array_of_slopes)
        print("######")
    return np.mean(array_of_slopes)

tfinal_list = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0]
RMS_list = []
RMS_list_transport = []
RMS_list_transport_2 = []
tfinal_list_3 = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0, 16384.0]
tfinal_list_2 = [16384.0]
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
def dipole_benchmark(xi):
    return 0.5 * (1/math.sqrt(math.pi)) * np.exp(-xi**2/4) * xi
def gaussian_IC_benchmark(xi):
    return 0.5  * np.exp(-xi**2/4)

def gaussian_IC_bench_full(z, tau, omega, A):
    exponential_term = -z**2/(omega**2 + 4* A*tau)
    return  np.exp(exponential_term) / np.sqrt(A * tau) / np.sqrt(4/omega**2 + 1/A/tau)




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



def ss_transport(N_spaces,scaled=True):
    RMS_list_transport = []
    sigma = 1
    order_list = []
    A = c/3/sigma
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
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize * sigma , benchmark(xi_sol)))
        if scaled == True:
            # if t == tfinal_list_2[0]:
            #     plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
            plt.plot(xi_sol, transport_phi / dimensionalize * sigma   , label = f't = {t}, {N_spaces} spatial cells')
        
        else:
            plt.plot(x, transport_phi, label = f't = {t}, {N_spaces} spatial cells')

        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)
        if scaled == True:
            plt.xlim(-10,10)
            plt.xlabel(r'$\xi$')
            plt.ylabel(r'$f(\xi)$')
            plt.legend(prop={'size':8})
        else:
            plt.xlabel('x')
            plt.ylabel(r'$\phi$')
            plt.xlim(-0.15, 0.15)
            plt.legend(prop={'size':8})



    # plt.legend()
    if scaled == True:
        show('self_similar_transport_pl')
    else:
        show('transport_not_scaled')
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    if t == tfinal_list_2[-1]:
        plt.loglog(tfinal_list_2, 0.07*np.array(tfinal_list_2)**-(np.sqrt(2)/2), label = r'$C_1\:t^{-\sqrt{2}}$')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.savefig('transport_convergence_pl.pdf')
    plt.show()
    

    if len(tfinal_list_2) > 1:
        print('### ### ###')
        print('average order', loglog_converge(np.array(tfinal_list_2), RMS_list_transport, 3), N_spaces, 'spatial cells')
        print('### ### ###')

    # return loglog_converge(np.array(tfinal_list_2), RMS_list_transport)
    return RMS_list_transport


def ss_special_IC(N_spaces):
    loader_transport = load_sol('transport', 'square_IC', 'transport', s2 = True, file_name = 'run_data.hdf5')
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
        # if t == tfinal_list_2[-1]:
        #     plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        plt.plot(xi_sol, transport_phi / dimensionalize * sigma   , label = f't = {t}, {N_spaces} spatial cells')
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize * sigma , benchmark(xi_sol)))
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        plt.xlim(-10,10)
        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$f(\xi)$')
        plt.legend(prop={'size':8})


    # plt.legend()
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.show()

def ss_dipole(N_spaces, scaled = True):
    RMS_list_transport = []
    sigma = 1600
    A = c/3/sigma
    for it, t in enumerate(tfinal_list_2):
        
        loader_transport.call_sol(t, 8, 0.25, N_spaces, 'rad', False, False)
        x0 = loader_transport.xs[-1]*2/N_spaces/2
        x = loader_transport.xs / sigma
        tau = t /c /sigma
        z0 = x0 / sigma
        t0 = -z0**2/6/A


        # t0 = 0
        # dimensionalize = x0*np.abs(x-x0) / ((A * (tau-t0))**1.5)
        dimensionalize = 1/(A * (tau-t0) /z0)
        l = x[-1] 
        early_time = 3 * sigma * z0**2 / c
        late_time = 0.1 * (l-z0)**2 / A

        early_time_nondim = early_time * c * sigma
        late_time_nondim = late_time * c * sigma
        # print(early_time_nondim, 'early time')
        if t == tfinal_list_2[-1]:
            print(late_time_nondim, 'late time')
            print(t0, 't0')
            print(z0, 'z0')
        
        transport_phi = loader_transport.phi 
        xi_sol = (x) / math.sqrt(A*(tau-t0))

        plt.figure(3)
        # plt.ion()
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        # if t == tfinal_list_2[-3]:
        #     plt.plot(xi_sol, dipole_benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize*sigma   , dipole_benchmark(xi_sol)))
        plt.ion()
        if scaled == True:
            plt.plot(xi_sol, transport_phi / dimensionalize* sigma     , label = f't = {t}')
            plt.xlim(0,5)
            plt.ylim(0,0.3)
        # plt.ylim(0,1)
            plt.xlabel(r'$\xi$')
            plt.ylabel(r'$f(\xi)$')
            plt.legend(prop={'size':8})
            show(f'dipole_scaled_{N_spaces}_cells')
        else:
            plt.plot(loader_transport.xs, transport_phi, label = f't = {t}')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$\phi$')
            plt.xlim(0, 30)
            plt.ylim(0, 0.075)
            plt.legend(prop={'size':8})
            show(f'dipole_notscaled_{N_spaces}_cells')
            


        # plt.plot(x*sigma, transport_phi     , label = f't = {t}, {N_spaces} spatial cells')
        
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        

        # plt.figure(4)
        # plt.plot(-x, transport_phi, label = f't = {t}, {N_spaces} spatial cells')
        # plt.xlim(0, 0.15)



    # plt.legend()
    plt.show()

    plt.figure(6)
    # plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells ' + r'$t_1=$' + '%0.0f' % late_time_nondim)
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.savefig('convergence_dipole')
    plt.show()
    print(loglog_converge(np.array(tfinal_list_2), RMS_list_transport), N_spaces, 'spatial cells')

def ss_transport_epsilon(N_spaces, epsilon = 1.0, scaled=True, tfinal_list_2 = [512.0]):
    RMS_list_transport = []
    sigma = 1
    order_list = []
    A = c/3/sigma
    for it, t in enumerate(tfinal_list_2):
        loader_transport.call_sol(t, 8, 1e-12, N_spaces, 'rad', False, False, epsilon)
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
        RMS_list_transport.append(RMSE(transport_phi / dimensionalize * sigma , benchmark(xi_sol)))
        if scaled == True:
            # if t == tfinal_list_2[0]:
            #     plt.plot(xi_sol, benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
            plt.plot(xi_sol, transport_phi / dimensionalize * sigma   , label = f't = {t}, {N_spaces} spatial cells')
        
        else:
            plt.plot(x, transport_phi, label = f't = {t}, {N_spaces} spatial cells')

        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize2 * sigma * 2 , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)
        if scaled == True:
            plt.xlim(-10,10)
            plt.xlabel(r'$\xi$')
            plt.ylabel(r'$f(\xi)$')
            plt.legend(prop={'size':8})
        else:
            plt.xlabel('x')
            plt.ylabel(r'$\phi$')
            plt.xlim(-0.15, 0.15)
            plt.legend(prop={'size':8})



    # plt.legend()
    if scaled == True:
        show('self_similar_transport_pl')
    else:
        show('transport_not_scaled')
    plt.show()

    plt.figure(6)
    plt.ion()
    plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells')
    # if t == tfinal_list_2[-1]:
    #     plt.loglog(tfinal_list_2, 0.07*np.array(tfinal_list_2)**-(np.sqrt(2)/2), label = r'$C_1\:t^{-\sqrt{2}}$')
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    # plt.savefig('transport_convergence_pl.pdf')
    plt.show()
    

    if len(tfinal_list_2) > 1:
        print('### ### ###')
        print('average order', loglog_converge(np.array(tfinal_list_2), RMS_list_transport, 3), N_spaces, 'spatial cells')
        print('### ### ###')

    # return loglog_converge(np.array(tfinal_list_2), RMS_list_transport)
    return RMSE(transport_phi / dimensionalize*sigma, gaussian_IC_benchmark(xi_sol))
def ss_gaussian_noniso(N_spaces, epsilon = 1.0, scaled = False, tfinal_list_2  = [1024.0]):
    loader_transport = load_sol('transport', 'gaussian_IC', 'transport', s2 = True, file_name = 'run_data.hdf5')
    RMS_list_transport = []
    sigma = 1
    A = c/3/sigma
    omega = 4.0
    for it, t in enumerate(tfinal_list_2):
    
        loader_transport.call_sol(t, 8, omega, N_spaces, 'rad', False, False, epsilon)
        x = loader_transport.xs / sigma
        tau = t /c /sigma
        transport_phi = loader_transport.phi
        plt.figure(3)
        dimensionalize = math.sqrt(1/omega**2) * math.sqrt(A*tau)
        xi_sol = x / math.sqrt(A*tau)

        # plt.ion()
        # plt.plot(xi_sol, loader.phi  * np.exp(-loader.xs**2/ 4/tau/c * 3* sigma)/ 2 * math.sqrt(math.pi * 3 * sigma/c*t) , label = f't = {t}')
        # plt.plot(xi, f, 'bo', mfc = 'none')
        # plt.plot(-xi, f, 'bo', mfc = 'none')
        # if t == tfinal_list_2[-3]:
        #     plt.plot(xi_sol, dipole_benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
        # RMS_list_transport.append())
        plt.ion()
        if scaled == True:
            plt.plot(xi_sol, transport_phi * dimensionalize *sigma -gaussian_IC_benchmark(xi_sol)  , label = r'$\epsilon =$'+ f'{epsilon}')
            # if t == tfinal_list_2[0]:
            #     plt.plot(xi_sol, gaussian_IC_benchmark(xi_sol), 'o', mfc = 'none', label = 'benchmark')
            # plt.xlim(0,5)
            # plt.ylim(0,0.3)
        # plt.ylim(0,1)
            plt.xlabel(r'$\xi$')
            plt.ylabel(r'$f(\xi)$')
            plt.legend(prop={'size':8})
            plt.xlim(-10,10)
            show(f'dipole_scaled_{N_spaces}_cells')
            return RMSE(transport_phi * dimensionalize*sigma, gaussian_IC_benchmark(xi_sol))
        else:
            plt.plot(x, transport_phi * sigma, label = f't = {t}')
            plt.plot(x, gaussian_IC_bench_full(x, tau, omega, A), 'o', mfc = 'none')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$\phi$')
            # plt.xlim(0, 30)
            # plt.ylim(0, 0.075)
            plt.legend(prop={'size':8})
            show(f'dipole_notscaled_{N_spaces}_cells')
    # 
    
            return RMSE(transport_phi * sigma, gaussian_IC_bench_full(x, tau, omega, A))
            


        # plt.plot(x*sigma, transport_phi     , label = f't = {t}, {N_spaces} spatial cells')
        
        # RMS_list_transport_2.append(RMSE(transport_phi / dimensionalize * sigma  , benchmark(xi_sol2)))
        # plt.plot(x, transport_phi)

        

        # plt.figure(4)
        # plt.plot(-x, transport_phi, label = f't = {t}, {N_spaces} spatial cells')
        # plt.xlim(0, 0.15)



    # plt.legend()
    plt.show()

    plt.figure(6)
    # plt.ion()
    # plt.loglog(tfinal_list_2, RMS_list_transport, '-o', label = f'{N_spaces} spatial cells ' + r'$t_1=$' + '%0.0f' % late_time_nondim)
    # plt.loglog(tfinal_list_2, RMS_list_transport_2, '-^', mfc = 'none')
    plt.xlabel('time')
    plt.ylabel('RMSE')
    plt.legend()
    plt.savefig('convergence_dipole')
    plt.show()
    # print(loglog_converge(np.array(tfinal_list_2), RMS_list_transport), N_spaces, 'spatial cells')


### transport plane ic convergence plot code block ###
# for it, tt in enumerate(tfinal_list_3):
#     tfinal_list_2 = [tt]
#     err_list = []
#     # order_list = []
#     N_spaces_list = [32, 64, 96, 128, 192, 256, 384, 512 ]
#     for NN in N_spaces_list:
#         err_list.append(ss_transport(NN))
#     plt.figure(123)
#     print('#################################################### ')
#     print(loglog_converge(np.array(N_spaces_list), np.array(err_list), 2), 'average order')
#     print('#################################################### ')
#     plt.loglog(N_spaces_list, err_list, '-o', mfc = 'none', label = (f't = {tt}'))
# # plt.loglog(N_spaces_list, 0.9*np.array(N_spaces_list)**(-np.sqrt(3)), label = r'$c_1\: K^{-\sqrt{3}}$')
# plt.loglog(N_spaces_list, 3.5*np.array(N_spaces_list)**(-2.), label = r'$c_1\: K^{-2}$')

# plt.rc('legend', fontsize=6)    # legend fontsize   
# plt.legend()     
# plt.xlabel(r'K')
# plt.ylabel('RMSE')
# plt.savefig('transport_converge_space.pdf')
# plt.close()
# plt.close()
# plt.close()
# plt.close()


### end block ###

# tfinal_list_2 = tfinal_list_3
# ss_transport(32)
# ss_transport(64)
# ss_transport(96)
# ss_transport(128)
# ss_transport(192)
# ss_transport(256)
# ss_transport(384)
# ss_transport(512)

# ss_suolson(32)
# ss_suolson(128)

# ss_dipole(32, True)
# ss_dipole(64, True)
# ss_dipole(96, True)
# ss_dipole(128, True)
# ss_dipole(192, True)
# ss_dipole(256)
# ss_dipole(384)


# ss_special_IC(64)
# ss_special_IC(32)

# ss_gaussian_noniso(16, 1.0)
# plt.close()


def epsilon_convergence_gauss(N_spaces = 32, epsilon_list = [0.5], scaled = True, tfinal_list = [1024.0]):
    rmse_list = []
    for iep, ep in enumerate(epsilon_list):
        rmse_list.append(ss_gaussian_noniso(N_spaces, ep, scaled = scaled, tfinal_list_2 = tfinal_list))
    
    plt.figure(11)
    plt.loglog(epsilon_list, rmse_list, '-o')
    shift = 0.494
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-0.002), label = r'$c_1\: \epsilon^{-0.002}$')
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-1), label = r'$c_1\: \epsilon^{-1}$')
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-3), label = r'$c_1\: \epsilon^{-3}$')
    plt.legend()
    plt.xlabel(r'$\epsilon$')
    plt.ylabel('RMSE')
    plt.show()
    print('###       convergence order            ###')
    print(loglog_converge(np.array(epsilon_list), rmse_list, 0))
    print('###                                    ###')

def epsilon_convergence_pl(N_spaces = 32, epsilon_list = [0.5], scaled = True, tfinal_list = [1024.0]):
    rmse_list = []
    for iep, ep in enumerate(epsilon_list):
        rmse_list.append(ss_transport_epsilon(N_spaces, ep, scaled = scaled, tfinal_list_2 = tfinal_list))
    
    plt.figure(11)
    plt.loglog(epsilon_list, rmse_list, '-o')
    shift = 0.494
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-0.002), label = r'$c_1\: \epsilon^{-0.002}$')
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-1), label = r'$c_1\: \epsilon^{-1}$')
    # plt.loglog(epsilon_list, shift*np.array(epsilon_list)**(-3), label = r'$c_1\: \epsilon^{-3}$')
    plt.legend()
    plt.xlabel(r'$\epsilon$')
    plt.ylabel('RMSE')
    plt.show()
    print('###       convergence order            ###')
    print(loglog_converge(np.array(epsilon_list), rmse_list, 0))
    print('###')              


