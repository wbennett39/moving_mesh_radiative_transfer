import numpy as np
from Chebyshev_matrix_reader import file_reader

class matrix_builder():
    def __init__(self, Mass_denom, J_denom, G_denom, L_denom, VV_denom, Mass_coeff_even, Mass_coeff_odd, 
                J_coeff_even, J_coeff_odd, G_coeff, L_coeff, VV_coeff_even, VV_coeff_odd):

        self.Mass_denom = Mass_denom[:]
        self.J_denom = J_denom[:]
        self.G_denom = G_denom[:]
        self.L_denom = L_denom[:]
        self.VV_denom = VV_denom[:]

        self.Mass_coeff_even = Mass_coeff_even[:]
        self.Mass_coeff_odd = Mass_coeff_odd[:]
        self.J_coeff_even = J_coeff_even[:]
        self.J_coeff_odd = J_coeff_odd[:]
        self.G_coeff = G_coeff
        self.L_coeff = L_coeff
        self.VV_coeff_even = VV_coeff_even
        self.VV_coeff_odd = VV_coeff_odd
        
        self.M = self.Mass_denom[0].size - 1
        self.N = self.VV_denom[:,0,0].size - 1
        self.Mass = np.zeros((self.M+1, self.M+1))
        self.J = np.zeros((self.M+1, self.M+1))
        self.G = np.zeros((self.M+1, self.M+1))
        self.L = np.zeros((self.M+1, self.M+1))
        self.VV = np.zeros((self.N+1, self.M+1, self.M+1))

    def make_mass(self, rL, rR):
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        pi = math.pi
        self.Mass[0,0] = (rL2 + rLrR + rR2) / 3 / pi
        # rL2 = rL**2
        # rR2 = rR**2
        # rL3 = rL**3
        # rR3 = rR**3
        # rLrR = rL*rR
        # even_fac = rL-rR
        # odd_fac = (rL-rR)**2

        # self.Mass[0,0] = self.Mass_coeff_odd[0, 0, 0] * rL3 + self.Mass_coeff_odd[0, 0, 1] * rR3 
        # for ii in range(self.M+1):
        #     for jj in range(self.M+1):
        #         if ii + jj != 0:
        #             if (ii + jj + 2) % 2 == 0 and ii + jj != 0:
        #                 self.Mass[ii][jj] = (self.Mass_coeff_even[ii, jj, 0] * rL2 + self.Mass_coeff_even[ii, jj, 1] * rLrR + 
        #                                     self.Mass_coeff_even[ii, jj, 2] * rR2) * even_fac
        #             else:
        #                 self.Mass[ii][jj] = (self.Mass_coeff_odd[ii, jj, 0] * rL + self.Mass_coeff_odd[ii, jj, 1] * rR )* odd_fac
        


        # self.Mass = np.multiply(self.Mass, 1/self.Mass_denom)

    def make_J(self, rL, rR):
        pi = math.pi
        self.J[0,0] = (rL + rR) / 2 /pi




        # rL2 = rL**2
        # rR2 = rR**2
        # rL3 = rL**3
        # rR3 = rR**3
        # rLrR = rL*rR
        # even_fac = rL-rR
        # odd_fac = (rL-rR)**2

        # self.J[0,0] = self.J_coeff_odd[0, 0, 0] * rL2 + self.J_coeff_odd[0, 0, 1] * rR2 
        # for ii in range(self.M+1):
        #     for jj in range(self.M+1):
        #         if ii + jj != 0:
        #             if (ii + jj + 2) % 2 == 0 and ii + jj != 0 :
        #                 self.J[ii][jj] = (self.J_coeff_even[ii, jj, 0] * rL + self.J_coeff_even[ii, jj, 1] * rR) * even_fac
        #             else:
        #                 self.J[ii][jj] =  odd_fac * self.J_coeff_odd[ii, jj, 0]
        


        # self.J = np.multiply(self.J, 1/self.J_denom)

    def make_G(self, rL, rR, rLp, rRp):
        rL2 = rL**2
        rR2 = rR**2
        rLrR = rL*rR
        pi = math.pi

        self.G[0,0] = - (rL2 + rLrR + rR2) * (rLp - rRp) / 6 / pi / (rL-rR)


        
        # rL2 = rL**2
        # rR2 = rR**2
        # rLrR = rL*rR

        # self.G[0,:] = 0.0
        # for ii in range(self.M+1):
        #     for jj in range(self.M+1):
        #         if ii > 0:
        #             self.G[ii,jj] = ((self.G_coeff[ii, jj, 0] * rL2 + self.G_coeff[ii, jj, 1] * rLrR + self.G_coeff[ii, jj, 2] * rR2) * rLp + 
        #                             (self.G_coeff[ii, jj, 3] * rL2 + self.G_coeff[ii, jj, 4] * rLrR + self.G_coeff[ii, jj, 5] * rR2) * rRp)


        # self.G[1:,:] = np.multiply(self.G[1:,:], 1/self.G_denom[1:,:])

    def make_L(self, rL, rR):
        self.L[0,0] = 0
        # rL2 = rL**2
        # rR2 = rR**2
        # rLrR = rL*rR

        # for ii in range(self.M+1):
        #     for jj in range(self.M+1):
        #         if ii > 0:
        #             self.L[ii,jj] = (self.L_coeff[ii, jj, 0] * rL2 + self.L_coeff[ii, jj, 1] * rLrR + self.L_coeff[ii, jj, 2] * rR2) 


        # self.L[1:,:] = np.multiply(self.L[1:,:], 1/self.L_denom[1:,:])
    
    def make_VV(self, rL, rR):
        # rL2 = rL**2
        # rR2 = rR**2
        # rL3 = rL**3
        # rR3 = rR**3
        # rLrR = rL*rR
        # even_fac = rL-rR
        # odd_fac = (rL-rR)**2

        # self.VV[0,0,0] = self.VV_coeff_odd[0, 0, 0,0] * rL3 + self.VV_coeff_odd[0,0, 0, 1] * rR3 
        # for nn in range(self.N+1):
        #     for ii in range(self.M+1):
        #         for jj in range(self.M+1):
        #             if nn + ii + jj != 0:
        #                 if (nn + 1 + ii + 1 + jj + 1 - 3) % 2 == 0 or (nn + jj + 2) %2 ==0:
        #                     self.VV[nn, ii, jj] = (self.VV_coeff_even[nn, ii, jj, 0] * rL2 + self.VV_coeff_even[nn, ii, jj, 1] * rLrR + 
        #                                         self.VV_coeff_even[nn, ii, jj, 2] * rR2) * even_fac
        #                 else:
        #                     self.VV[nn, ii, jj] = (self.VV_coeff_odd[nn, ii, jj, 0] * rL + self.VV_coeff_odd[nn, ii, jj, 1] * rR)* odd_fac
        

        # self.VV = np.multiply(self.VV, 1/self.VV_denom)
        self.VV[0,0] = 0
        
        
        

def mass_builder_test():
    # works for M = 3
    reader = file_reader()
    give = reader()
    ob = matrix_builder(give[0], give[1], give[2], give[3], give[4], give[5], give[6], give[7], give[8], give[9], give[10], give[11], give[12])
    ob.make_mass(1.246, 3.2349)
    # print(1/ob.Mass_denom)
    if ob.M == 3:
        print('testing mass matrix')
        answer_matrix = np.array([[10.6391, 1.31281, -2.82287, -0.737613], [1.31281, 
                                1.71941, 0.287596, -1.00307], [-2.82287, 0.287596, 
                                5.06152, 0.556503], [-0.737613, -1.00307, 0.556503, 
                                2.45695]])
        np.testing.assert_allclose(ob.Mass, answer_matrix, rtol = 1e-5)


        answer_matrix = np.array([[5334.34, 1333.58, -1981.32, -800.15], [1333.58, 1168.47, 
                                    266.717, -593.651], [-1981.32, 266.717, 2526.3, 571.536], 
                                    [-800.15, -593.651, 571.536, 1538.95]])
        ob.make_mass(0, 25.2)
        np.testing.assert_allclose(ob.Mass, answer_matrix, rtol = 1e-5)
        

def J_builder_test():
    # works for M = 3
    reader = file_reader()
    give = reader()
    ob = matrix_builder(give[0], give[1], give[2], give[3], give[4], give[5], give[6], give[7], give[8], give[9], give[10], give[11], give[12])
    ob.make_J(0.0112, 2.3)
    # print(1/ob.Mass_denom)
    if ob.M == 3:
        print('testing J matrix')
        answer_matrix = np.array([[2.64494, 0.34924, -0.87481, -0.182935], [0.34924, 
                                0.439236, 0.0831525, -0.263743], [-0.87481, 0.0831525, 
                                1.23496, 0.146418], [-0.182935, -0.263743, 0.146418, 
                                0.640295]])
        np.testing.assert_allclose(ob.J, answer_matrix, rtol = 1e-5)
        
def G_builder_test():
    # works for M = 3
    reader = file_reader()
    give = reader()
    ob = matrix_builder(give[0], give[1], give[2], give[3], give[4], give[5], give[6], give[7], give[8], give[9], give[10], give[11], give[12])
    if ob.M == 3:
        print('Testing G matrix')
        ob.make_G(0.42, 1.2, 0.1, 2.3)
        answer_matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0], [-2.77708, -1.36045, 
                        0.61974, 0.799098, 0.356574, 0.204476, 0.12592, 
                        0.093148], [-4.32251, -3.81712, -0.868101, 1.88638, 
                        1.59962, 0.765048, 0.451086, 0.313371], [-3.61769, 
                        -4.44236, -4.97099, -1.24353, 3.04051, 2.76449, 1.37603, 
                        0.872109], [-3.4724, -3.8615, -5.44577, -6.10415, 
                        -0.834028, 4.39949, 3.65833, 1.88573], [-3.39215, 
                        -3.61389, -4.249, -6.35325, -7.45529, -1.23442, 5.57079, 
                        4.87138], [-3.36978, -3.4971, -3.85534, -4.85213, 
                        -7.48003, -8.62277, -0.828178, 6.94605], [-3.35021, 
                        -3.44807, -3.67671, -4.26331, -5.35802, -8.45159, 
                        -9.99588, -1.23205]])

       
        np.testing.assert_allclose(ob.G, answer_matrix[0:4,0:4], rtol = 1e-5)
        

def L_builder_test():
    # works for M = 3
    reader = file_reader()
    give = reader()
    ob = matrix_builder(give[0], give[1], give[2], give[3], give[4], give[5], give[6], give[7], give[8], give[9], give[10], give[11], give[12])
    if ob.M == 3:
        print('Testing L matrix')
        ob.make_L(4.2, 6.2)
        answer_matrix =  np.array([[0, 0, 0, 0], [37.9067, 18.0067, -18.328, -11.029], 
                                [23.2533, 56.641, 5.33333, -35.7705], [38.7194, 51.0568, \
                                73.9869, 17.365]])

       
        np.testing.assert_allclose(ob.L, answer_matrix[0:4,0:4], rtol = 1e-5)
def VV_builder_test():
    # works for M = 3
    reader = file_reader()
    give = reader()
    ob = matrix_builder(give[0], give[1], give[2], give[3], give[4], give[5], give[6], give[7], give[8], give[9], give[10], give[11], give[12])
    if ob.M == 3:
        print('Testing L matrix')
        ob.make_VV(4.2, 6.2)
        answer_matrix =  np.array([[[54.7467, 2.90667, -8.94248, -1.57333], [2.90667, 7.08013, 
                                    0.666667, -4.47132], [-8.94248, 0.666667, 26.4629, 
                                    1.22483], [-1.57333, -4.47132, 1.22483, 10.6005]], 
                                    [[2.90667, 7.08013, 0.666667, -4.47132], [7.08013, 
                                    3.88036, 1.30441, -0.573945], [0.666667, 1.30441, 
                                    0.945748, 3.06459], [-4.47132, -0.573945, 3.06459, 
                                    2.91088]], [[-8.94248, 0.666667, 26.4629, 1.22483], 
                                    [0.666667, 1.30441, 0.945748, 3.06459], [26.4629, 
                                    0.945748, -6.90314, -0.11557], [1.22483, 3.06459, 
                                    -0.11557, -4.2126]], [[-1.57333, -4.47132, 1.22483, 
                                    10.6005], [-4.47132, -0.573945, 3.0645, 2.91088], 
                                    [1.22483, 3.06459, -0.11557, -4.2126], [10.6005, 2.91088, 
                                    -4.2126, -3.01607]], [[-1.82084, -1.01517, -4.86381, 
                                    1.34219], [-1.01517, -2.71113, 0.163511, 3.32331], 
                                    [-4.86381, 0.163511, 12.667, 0.186086], [1.34219, 
                                    3.32331, 0.186086, -2.62021]]])

       

        res = np.abs(answer_matrix-ob.VV)
        for nn in range(5):
            for ii in range(4):
                for jj in range(4):
                    if res[nn, ii, jj] >= 1e-5:
                        print(nn, ii, jj)
                        print(answer_matrix[nn, ii, jj])
                        print(ob.VV[nn, ii, jj])
        np.testing.assert_allclose(ob.VV, answer_matrix, rtol = 1e-4)

                    

                   
# mass_builder_test()
# J_builder_test()
# G_builder_test()
# L_builder_test()
# VV_builder_test()