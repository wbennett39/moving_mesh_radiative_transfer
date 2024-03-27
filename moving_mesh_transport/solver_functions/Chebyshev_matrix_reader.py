import h5py 
import numpy as np

class file_reader():
    def __init__(self, mass_matrix_name = 'moving_mesh_transport/matrix_files/MassMatrix.h5', J_matrix_name = 'moving_mesh_transport/matrix_files/JMatrix.h5', G_matrix_name = 'moving_mesh_transport/matrix_files/GMatrix.h5',
        L_matrix_name = 'moving_mesh_transport/matrix_files/LMatrix.h5', VV_matrix_name = 'moving_mesh_transport/matrix_files/VVMatrix.h5'):
        """ Read in the mass matrix coefficients and denominators
        """
        f = h5py.File(mass_matrix_name, 'r+')
        f1 = h5py.File(J_matrix_name, 'r+')
        f2 = h5py.File(G_matrix_name, 'r+')
        f3 = h5py.File(L_matrix_name, 'r+')
        f4  = h5py.File(VV_matrix_name, 'r+')
        self.mass_mat_denom = f['denominators'][:]
        self.J_mat_denom = f1['denominators'][:]
        self.G_mat_denom = f2['denominators'][:]
        self.L_mat_denom = f3['denominators'][:]
        self.VV_mat_denom = f4['denominators'][:]
        M = self.mass_mat_denom[0].size
        # M =2
        N = self.VV_mat_denom[:,0,0].size
        self.mass_mat_coeff_even = np.zeros((M, M, 3))
        self.mass_mat_coeff_odd = np.zeros((M, M, 2))
        self.J_mat_coeff = np.zeros((M, M, 2))
        self.G_mat_coeff = np.zeros((M, M, 6))
        self.L_mat_coeff_even = np.zeros((M, M, 2))
        self.L_mat_coeff_odd = np.zeros((M, M, 3))
        self.VV_mat_coeff_even = np.zeros((N, M, M, 3))
        self.VV_mat_coeff_odd = np.zeros((N, M, M, 2))
        for ii in range(M):
            for jj in range(M):
                self.J_mat_coeff[ii, jj] = f1['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
                # if ii > 0:
                # self.G_mat_coeff[ii, jj] = f2['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
                
                if (ii + jj + 2) % 2 == 0:
                    self.mass_mat_coeff_even[ii, jj] = f['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
                    if ii > 0:
                        self.L_mat_coeff_even[ii, jj] = f3['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
                    
                else:
                    self.mass_mat_coeff_odd[ii, jj] = f['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
                    if ii >0:
                        self.L_mat_coeff_odd[ii, jj] = f3['coefficients'][f'Element{ii+1}{jj+1}'][:].flatten()
        # for nn in range(N):
        #     for ii in range(M):
        #         for jj in range(M):
        #             if ((nn + 1 + ii + 1 + jj + 1 - 3) % 2 == 0 and nn + ii + jj != 0) or ((nn + jj + 2) % 2 == 0  and nn + ii + jj != 0):
        #                 self.VV_mat_coeff_even[nn, ii, jj] = f4['coefficients'][f'Element{nn+1}{ii+1}{jj+1}'][:].flatten()
        #             else:
        #                 self.VV_mat_coeff_odd[nn, ii, jj] = f4['coefficients'][f'Element{nn+1}{ii+1}{jj+1}'][:].flatten()

        f.close()
        f1.close()
        f2.close()
        f3.close()
        f4.close()

    def __call__(self):
        return  self.mass_mat_denom, self.J_mat_denom, self.G_mat_denom, self.L_mat_denom, self.VV_mat_denom, self.mass_mat_coeff_even, self.mass_mat_coeff_odd, self.J_mat_coeff, self.G_mat_coeff, self.L_mat_coeff_even, self.L_mat_coeff_odd, self.VV_mat_coeff_even, self.VV_mat_coeff_odd



            
