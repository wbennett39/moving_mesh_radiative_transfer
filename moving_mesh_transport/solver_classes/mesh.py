import numpy as np
from numba import int64, float64
import numba
from numba.experimental import jitclass
import math
from .functions import problem_identifier 
from .mesh_functions import set_func, _interp1d
import quadpy
import numpy.polynomial as nply
from scipy.special import roots_legendre


#################################################################################################
data = [('N_ang', int64), 
        ('N_space', int64),
        ('M', int64),
        ('tfinal', float64),
        ('mus', float64[:]),
        ('ws', float64[:]),
        ('x0', float64),
        ("moving", int64),
        ("move_type", int64[:]),
        ("edges", float64[:]),
        ("edges0", float64[:]),
        ("Dedges", float64[:]),
        ("N_space", int64),
        ('middlebin', int64),
        ('sidebin', int64),
        ('speed', float64),
        ('Dedges_const', float64[:]),
        ('source_type', int64[:]),
        ('thick', int64), 
        ('move_func', int64),
        ('debugging', int64),
        ('wave_loc_array', float64[:,:,:]),
        ('delta_t', float64),
        ('tactual', float64),
        ('told', float64),
        ('index_old', int64),
        ('right_speed', float64),
        ('left_speed', float64),
        ('test_dimensional_rhs', int64),
        ('move_factor', float64),
        ('T_wave_speed', float64),
        ('pad', float64),
        ('follower_speed', float64),
        ('leader_speed', float64),
        ('span_speed', float64),
        ('thick_quad', float64[:]),
        ('middlebin', int64),
        ('sidebin', int64),
        ('leader_pad', float64),
        ('packet_leader_speed', float64),
        ('thick_quad_edge', float64[:]),
        ('t0', float64),
        ('edges0_2', float64[:])
        # ('problem_type', int64)
        ]
#################################################################################################


    
@jitclass(data)
class mesh_class(object):
    def __init__(self, N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick, move_factor, wave_loc_array, pad, leader_pad, thick_quad, thick_quad_edge):
        
        self.debugging = True
        self.test_dimensional_rhs = False

        self.pad = pad
        self.tfinal = tfinal
        self.N_space = N_space
        self.x0 = x0
        self.moving = moving
        self.move_type = move_type
        self.edges = np.zeros(N_space+1)
        self.edges0 = np.zeros(N_space+1)
        self.Dedges = np.zeros(N_space+1)
        self.N_space = N_space
        self.speed = edge_v
        self.move_factor = move_factor
        if self.test_dimensional_rhs == True:
            self.speed = 299.98

        print('mesh edge velocity: ', edge_v)
        self.source_type = source_type

        if self.move_type[0] == True:
            self.move_func = 0 # simple linear
        elif self.move_type[1] == True:
            self.move_func = 1 # thick square source move
            print('thick square source edge estimation mesh')
        elif self.move_type[2] == True:
            self.move_func = 2 # sqrt t static
        elif self.move_type[3] == True:
            self.move_func = 3 # sqrt t static
        
        # self.problem_type = problem_identifier(self.source_typem, self.x0)
        self.thick = thick
        self.wave_loc_array = wave_loc_array
        # self.smooth_wave_loc_array
        # for count, element in enumerate(self.wave_loc_array[0, 3, :]):
        #     if element < self.x0:
        #         self.wave_loc_array[0, 3, count] = self.x0 + 1e-8
        self.thick_quad = thick_quad
        self.thick_quad_edge = thick_quad_edge

        # print(self.wave_loc_array[0,2,-1], 'wave final location')

        self.sidebin = int(self.N_space/4)
        self.middlebin = int(self.N_space/2)
        
        self.initialize_mesh()
        self.tactual = -1.
        self.told = 0.0
        self.index_old = 0
        self.T_wave_speed = 0.0
        self.follower_speed = 0.0
        self.leader_speed = 0.0
        self.span_speed = 0.0
        self.leader_pad = leader_pad
        self.t0 = 10.0
    

    def move(self, t):
        # print(self.edges)pr
        """
        Called each time the rhs moves the mesh. Changes edges and Dedges
        """
        # print(self.edges)
        # if self.moving == True:
        """
        This mode moves all of the edges at a constant speed,
        linearly increasing from 0 to the wavespeed
        """
        if self.moving == True:
            # if self.source_type[1] == 1 or self.source_type[2] == 1:
                # if t > 10.0:
                #     self.Dedges = self.edges/self.edges[-1] * self.speed
            if self.move_func == 0: # simple linear
                if t > self.t0 and self.source_type[2] == 1:
                    self.move_middle_edges(t)
                    self.edges = self.edges0_2 + self.Dedges * (t-self.t0)
     

                    # print(self.Dedges)


                if t <= self.t0 or self.source_type[2] != 1:
                    self.Dedges = self.Dedges_const
                    self.edges = self.edges0 + self.Dedges*t
                else:
                    self.edges = self.edges0 + self.Dedges*t

            elif self.move_func == 1: 
                """
                This mode has the wavefront tracking edges moving at a constant speed
                and interior edges tracking the diffusive wave
                """
                # self.thick_square_moving_func(t)
                self.thick_square_moving_func_2(t)
                

            elif self.move_func == 2:
                self.square_source_static_func_sqrt_t(t)


        

            else:
                print("no move function selected")
                assert(0)

            # if self.debugging == True:
            #     for itest in range(self.edges.size()):
            #         if self.edges[itest] != np.sort(self.edges)[itest]:
            #             print("crossed edges")
            #             assert(0)
    def smooth_wave_loc_array(self):
        for ix in range(0,self.wave_loc_array[0,3,:].size-1):
            if self.wave_loc_array[0,3,ix] < self.wave_loc_array[0,3,ix +1]:
                self.wave_loc_array[0,3,ix] = self.wave_loc_array[0,3,ix +1]

    def move_middle_edges(self,t):
        middlebin = int(self.N_space/2)
        sidebin = int(middlebin/2)
        if self.Dedges[sidebin] == 0:
            # final_pos = self.edges0[-1] + self.Dedges[-1] * self.tfinal
            final_pos = self.pad
            # final_pos = self.x0 + self.pad
            final_array = np.linspace(-final_pos, final_pos, self.N_space + 1)
            new_Dedges = (final_array - self.edges) / (self.tfinal-self.t0)
            # self.Dedges = new_Dedges
            self.Dedges = self.Dedges_const
            self.edges0_2 = self.edges
            
            # loc_of_right_outside_edge = self.edges[middlebin+sidebin] + self.Dedges[sidebin] * (self.tfinal-t)
            # loc_of_left_outside_edge = self.edges[sidebin] + self.Dedges[sidebin] * (self.tfinal-t)
            # dx = loc_of_right_outside_edge - loc_of_left_outside_edge
            # speed_right = dx/2/(self.tfinal-t)
            # speed_array = self.edges[sidebin:middlebin+sidebin+1]/speed_right
            # self.Dedges[sidebin:sidebin+middlebin+1] = speed_array
            # print(' --- --- --- --- --- --- --- ')
            # print('speed_array')
            # print(speed_array)
            # print(' --- --- --- --- --- --- --- ')

    

    def thick_wave_loc_and_deriv_finder(self, t):
        
        interpolated_wave_locations = _interp1d(np.ones(self.wave_loc_array[0,3,:].size)*t, self.wave_loc_array[0,0,:], self.wave_loc_array[0,3,:], np.zeros(self.wave_loc_array[0,3,:].size))
        # interpolated_wave_locations = np.interp(t, self.wave_loc_array[0,0,:], self.wave_loc_array[0,3,:] )

        # derivative = (interpolated_wave_locations[1] - interpolated_wave_locations[0])/delta_t_2

        edges = np.copy(self.edges)
        if t == 0 or  interpolated_wave_locations[0] < self.x0:
            edges = self.edges0
        else:
            edges[-1] = interpolated_wave_locations[0] + self.leader_pad
            if edges[-1] < self.edges0[-1]:
                edges[-1] = self.edges0[-1]
            edges[-2] = interpolated_wave_locations[0] + self.pad
            if edges[-2] < self.edges0[-2]:
                edges[-2] = self.edges0[-2]
            edges[-3] = interpolated_wave_locations[0] 
            if edges[-3] < self.edges0[-3]:
                edges[-3] = self.edges0[-3]
            edges[-4] = interpolated_wave_locations[0]  - self.pad
            if edges[-4] < self.edges0[-4]:
                edges[-4] = self.edges0[-4]

        
            edges[self.sidebin + self.middlebin: self.N_space-3] = np.linspace(self.x0, edges[-4], self.sidebin - 2)[:-1] 
            edges[0:self.sidebin] =  - np.flip(np.copy(edges[self.sidebin+self.middlebin+1:]))
        return edges

    def thick_square_moving_func_2(self, t):
        delta_t = 1e-7
        self.edges = self.thick_wave_loc_and_deriv_finder(t)
        edges_new = self.thick_wave_loc_and_deriv_finder(t + delta_t)

        self.Dedges = (edges_new - self.edges) / delta_t

        if self.edges[-3] < self.edges[-4]:
            print('crossed')

    def recalculate_wavespeed(self, t):
        sidebin = int(self.N_space/4)
        T_index = -2
        index = np.searchsorted(self.wave_loc_array[0,0,:], t)
        # pad = self.edges[int(self.N_space/2)+1] - self.edges[int(self.N_space/2)
        # print(abs(self.edges[-2]-self.wave_loc_array[0,3,index+1]), 'T wave from mesh edge')
        if self.debugging == True:
            if index >0:
                if not (self.wave_loc_array[0,0,index-1] <= t <= self.wave_loc_array[0,0,index+1]):
                    print('indexing error')
        if index != self.index_old:
            # print('calculating wavespeeds') s
            T_wave_location = self.wave_loc_array[0,3,index+1]
            # print(self.pad, 'pad')
            # print(index, 'index')
            # print(T_wave_location, 'T wave loc')
            self.delta_t = self.wave_loc_array[0,0,index+1] - t
            self.right_speed = (self.wave_loc_array[0,2,index+1]  - self.edges[-1])/self.delta_t
            self.T_wave_speed = (T_wave_location - self.edges[-2])/self.delta_t
            # print(T_wave_location, 't edge is moving to')
            self.leader_speed = (T_wave_location + self.leader_pad - self.edges[-1])/self.delta_t
            self.packet_leader_speed = (T_wave_location + self.pad - self.edges[-2])/self.delta_t
            # print(T_wave_location + self.pad, 'leader edge is moving to')
            self.follower_speed = (T_wave_location - self.pad - self.edges[-3])/self.delta_t

            last_follower_edge_loc = self.edges[-3] + self.Dedges_const[-3] * self.follower_speed * self.delta_t
            dx_span = (last_follower_edge_loc - self.x0) / (sidebin/2)  
            self.span_speed = (last_follower_edge_loc - dx_span - self.edges[-int(sidebin-2)])/self.delta_t
        
        self.index_old = index
        # print(self.edges)

        # print(self.delta_t, 'delta t')
        # print(self.leader_speed, 'leader')
        # print(self.T_wave_speed, "T speed")
        # print(self.follower_speed, 'follower')
        # # print(self.span_speed, 'span')
        # print(self.T_wave_speed, 't wave s')
        # print(self.leader_speed, 'leader s')
        # print(index, 'index')
        if self.T_wave_speed > self.leader_speed:
            print("speed problem")
            print(self.pad, 'pad')

      
    
        # if abs(self.edges[T_index] + self.Dedges_const[T_index] * self.T_wave_speed * self.delta_t - (self.edges[T_index -1] + self.Dedges_const[T_index-1] * self.T_wave_speed * self.delta_t)) <= 1e-12:
        #     #catching up
        #     print('catching up')
        #     self.T_wave_speed = (self.edges[(T_index)-1] - 0.0005 - self.edges[(T_index)])/self.delta_t

        
        if self.debugging == True:
            if abs(t - self.wave_loc_array[0,0,index+1]) < 1e-5:
                print('checking location')
                print(self.wave_loc_array[0,3,index+1] - self.edges[-2], 'T wave difference')
                print(self.wave_loc_array[0,2,index+1]  - self.edges[-1], 'right edge difference')
        

        if self.right_speed < 0.0:
            self.right_speed = 0.0
        if self.T_wave_speed < 0.0:
            # print('negative t speed')
            self.T_wave_speed = 0.0
        if self.follower_speed < 0.0:
            self.follower_speed = 0.0
        if self.leader_speed < 0.0:
            self.leader_speed = 0.0
        if self.span_speed < 0.0:
            self.span_speed = 0.0
        
        
        # print(self.edges[-2], "|", self.wave_loc_array[0,3,index+1])


    def initialize_mesh(self):
        """
        Initializes initial mesh edges and initial edge derivatives. This function determines
        how the mesh will move
        """
        # if self.problem_type in ['plane_IC']:
        if self.source_type[0] == 1:
            self.simple_moving_init_func()

        if self.thick == False:     # thick and thin sources have different moving functions

            # if self.problem_type in ['gaussian_IC', 'gaussian_source']:
            if self.source_type[3] == 1 or self.source_type[5] == 1:
                self.simple_moving_init_func()
            # elif self.problem_type in ['square_IC', 'square_source']:
            if self.source_type[1] == 1 or self.source_type[2] == 1:
                self.thin_square_init_func_legendre()

        elif self.thick == True:
            # if self.problem_type in ['gaussian_IC', 'gaussian_source']:
            if self.source_type[3] == 1 or self.source_type[5] == 1:
                if self.moving == True:
                    self.simple_moving_init_func()
                elif self.moving == False:
                    self.thick_gaussian_static_init_func()

            elif self.source_type[1] == 1 or self.source_type[2] == 1:
                if self.move_func == 0:
                    self.simple_moving_init_func()

                if self.moving == False:
                    if self.move_func == 1:
                        self.simple_thick_square_init_func()
                    elif self.move_func == 2:
                        self.thin_square_init_func()
                elif self.moving == True:
                    if self.move_func == 1:
                        self.thick_square_moving_init_func()
                    elif self.move_func == 2:
                        self.thin_square_init_func()

        self.edges0 = np.copy(self.edges)
        self.Dedges_const = np.copy(self.Dedges)



        if self.moving == False:
            self.tactual = 0.0
            # static mesh -- puts the edges at the final positions that the moving mesh would occupy
            # sets derivatives to 0
            self.moving = True
            if self.thick == True:
                self.delta_t = self.tfinal 
            self.move(self.tfinal)
            self.Dedges = self.Dedges*0
            self.moving = False
            # self.edges[-1] = self.x0 + self.tfinal * self.speed
            # self.edges[0] = -self.x0 + -self.tfinal * self.speed


            print(self.edges[-1], "final edges -- last edge")

            

    def thick_gaussian_static_init_func(self):
        # if abs(self.wave_loc_array[0, 2, -1]) > 5:
        if self.move_func == 1:
            right_edge = self.wave_loc_array[0,3,-1] + self.pad
        elif self.move_func == 0:
            right_edge = self.x0
        print(self.move_func, 'move_func')
        print(right_edge, 'right edge')
        # else:
            # right_edge = self.x0 + self.tfinal
        
        # if right_edge < self.x0:
            # right_edge = self.x0 + self.tfinal

        self.edges = np.linspace(-right_edge, right_edge, self.N_space + 1)
        self.Dedges = self.edges * 0


    def simple_thick_square_init_func(self):
        # does not accomodate moving mesh edges
        
        # wave_edge = self.wave_loc_array[0,2,index+1]
        wave_edge = self.wave_loc_array[0,2,-1] + self.pad

        if self.N_space == 2:
            print("don't run this problem with 2 spaces")
            assert(0)
        middlebin = int(self.N_space/2)   # edges inside the source - static
        sidebin = int(middlebin/2) # edges outside the source - moving
        dx = 1e-8
        left = np.linspace(-wave_edge, -self.x0, sidebin + 1)
        right = np.linspace(self.x0, wave_edge, sidebin + 1)
        middle = np.linspace(-self.x0, self.x0, middlebin + 1)
        self.edges = np.concatenate((left[:-1], middle[:-1], right[:])) # put them all together 
        
        # initialize derivatives
        self.Dedges[0:sidebin] = (self.edges[0:sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        self.Dedges[sidebin:sidebin+middlebin] = 0       
        self.Dedges[middlebin+sidebin + 1:] = (self.edges[middlebin+sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
        self.Dedges = self.Dedges * self.speed * 0 

        
    def square_source_static_func_sqrt_t(self, t):
        # only to be used to estimate the wavespeed
        move_factor = self.move_factor

        if t > 1e-10:
            sqrt_t = math.sqrt(t)
        else:
            sqrt_t = math.sqrt(1e-10)

        # move the interior edges
        self.Dedges = self.Dedges_const * move_factor * 0.5 / sqrt_t
        self.edges = self.edges0 + self.Dedges_const * move_factor * sqrt_t
    

        # move the wavefront edges
        # Commented code below moves the exterior edges at constant speed. Problematic because other edges pass them
        # self.Dedges[0] = self.Dedges_const[0]
        # self.Dedges[-1] = self.Dedges_const[-1]
        # self.edges[0] = self.edges0[0] + self.Dedges[0]*t
        # self.edges[-1] = self.edges0[-1] + self.Dedges[-1]*t
        # self.Dedges[0] = self.Dedges_const[0] * move_factor * 0.5 / sqrt_t
        # self.edges[0] = self.edges0[0] + self.Dedges_const[0] * move_factor * sqrt_t
        # print(self.edges0[0], 'x0')
        # print(self.Dedges_const[0]*move_factor*sqrt_t, 'f')
        # self.Dedges[-1] = self.Dedges_const[-1] * move_factor * 0.5 / sqrt_t
        # self.edges[-1] = self.edges0[-1] + self.Dedges_const[-1] * move_factor * sqrt_t

    ####### Initialization functions ########


    def simple_moving_init_func(self):
            self.edges = np.linspace(-self.x0, self.x0, self.N_space+1)
            self.Dedges = self.edges/self.edges[-1] * self.speed

    # def thick_square_moving_func(self, t):
    #     middlebin = int(self.N_space/2)
    #     sidebin = int(self.N_space/4)
    #     self.recalculate_wavespeed(t)
    #     little_delta_t = t-self.told

    #     # self.Dedges[0:sidebin/2] = self.Dedges_const[0:sidebin/2] * self.right_speed
    #     # self.Dedges[0:int(sidebin/4 + 1)] = self.Dedges_const[0:int(sidebin/4 + 1)] * self.leader_speed
    #     # self.Dedges[int(sidebin/4 + 1)] = self.Dedges_const[int(sidebin/4 + 1)] * self.T_wave_speed
    #     # self.Dedges[int(sidebin/4 + 2):int(sidebin/2 + 1)] = self.Dedges_const[int(sidebin/4 + 2):int(sidebin/2 + 1)] * self.follower_speed
    #     # self.Dedges[int(sidebin/2 + 2):sidebin] = self.Dedges_const[int(sidebin/2 + 2):sidebin] * self.span_speed

    #     self.Dedges[0] = self.Dedges_const[0] * self.leader_speed
    #     self.Dedges[1] = self.Dedges_const[1] * self.T_wave_speed
    #     self.Dedges[2] = self.Dedges_const[2] * self.follower_speed
    #     self.Dedges[3:sidebin] = self.Dedges_const[3:sidebin] * self.span_speed

    #     self.Dedges[middlebin+sidebin + 1:] =  - np.flip(np.copy(self.Dedges[0:sidebin]))

    #     # self.Dedges[1:-1] =  self.Dedges_const[1:-1] * self.T_wave_speed 
    #     # self.Dedges[0] =  self.Dedges_const[0] * self.right_speed 
    #     # self.Dedges[-1] =  self.Dedges_const[-1] * self.right_speed 


    #     # self.Dedges = self.Dedges_const * self.right_speed
    #     # self.edges = self.edges + self.Dedges * delta_t
    #     # self.told = t
    #     # print(self.edges[-1]-self.edges[-2], 'thin zone')

    #     self.edges = self.edges + self.Dedges * little_delta_t
    #     self.told = t

    def thin_square_init_func(self):
        if self.N_space == 2:
            print("don't run this problem with 2 spaces")
            assert(0)
        middlebin = int(self.N_space/2)   # edges inside the source - static
        sidebin = int(middlebin/2) # edges outside the source - moving
        dx = 1e-14
        left = np.linspace(-self.x0-dx, -self.x0, sidebin + 1)
        right = np.linspace(self.x0, self.x0 + dx, sidebin + 1)

        
        middle = np.linspace(-self.x0, self.x0, middlebin + 1)

        self.edges = np.concatenate((left[:-1], middle[:-1], right[:])) # put them all together 
        
        # initialize derivatives
        self.Dedges[0:sidebin] = (self.edges[0:sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        self.Dedges[sidebin:sidebin+middlebin] = 0       
        self.Dedges[middlebin+sidebin + 1:] = (self.edges[middlebin+sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
        self.Dedges = self.Dedges * self.speed


    def thin_square_init_func_legendre(self):
        if self.N_space == 2:
            print("don't run this problem with 2 spaces")
            assert(0)
        middlebin = int(self.N_space/2)   # edges inside the source - static
        sidebin = int(middlebin/2) # edges outside the source - moving
        dx = 1e-8
        # left = np.linspace(-self.x0-dx, -self.x0, sidebin + 1)
        # right = np.linspace(self.x0, self.x0 + dx, sidebin + 1)
        left_old = self.thick_quad_edge
        right_old = self.thick_quad_edge
        right =(right_old*(self.x0-self.x0-dx)-self.x0-dx-self.x0)/-2
        left =(left_old*(-self.x0-dx+self.x0)+self.x0+dx+self.x0)/-2

        # if self.N_space == 32 and self.move_func == 2:
        #     middle = np.array([-0.99057548, -0.95067552, -0.88023915, -0.781514  , -0.65767116,
        #                 -0.51269054, -0.35123176, -0.17848418,  0.        ,  0.17848418,
        #                 0.35123176,  0.51269054,  0.65767116,  0.781514  ,  0.88023915,
        #                 0.95067552,  0.99057548]) 
        # else:
        # if self.move_func == 2:
        middle = 0.5 * self.thick_quad

            # left = roots_legendre(siebin+1)[0]
            # right = roots_legendre(siebin+1)[0]
            # right =(right*(self.x0-self.x0-dx)-self.x0-dx-self.x0)/-2
            # left =(left*(-self.x0-dx+self.x0)+self.x0+dx+self.x0)/-2
            # print(left, right)

        self.edges = np.concatenate((left[:-1], middle[:-1], right[:])) # put them all together 
        
        # initialize derivatives
        self.Dedges[0:sidebin] = (self.edges[0:sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        self.Dedges[sidebin:sidebin+middlebin] = 0       
        self.Dedges[middlebin+sidebin + 1:] = (self.edges[middlebin+sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
        self.Dedges = self.Dedges * self.speed


    def simple_thick_square_init_func_2(self):
        if self.N_space == 2:
            print("don't run this problem with 2 spaces")
            assert(0)
        middlebin = int(self.N_space/2)   # edges inside the source - static
        sidebin = int(middlebin/2) # edges outside the source - moving
        dx = 1e-14
        left = np.linspace(-self.x0-dx, -self.x0, sidebin + 1)
        right = np.linspace(self.x0, self.x0 + dx, sidebin + 1)
        middle = np.linspace(-self.x0, self.x0, middlebin + 1)
        self.edges = np.concatenate((left[:-1], middle[:-1], right[:])) # put them all together 
        
        # initialize derivatives
        self.Dedges[0:sidebin] = (self.edges[0:sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        self.Dedges[sidebin:sidebin+middlebin] = 0       
        self.Dedges[middlebin+sidebin + 1:] = (self.edges[middlebin+sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)
        self.Dedges = self.Dedges * self.speed * 0
    
    def thick_square_moving_init_func(self):
        # if self.N_space ==2 or self.N_space == 4 or self.N_space == 8:
        #     print(f"don't run this problem with {self.N_space} spaces")
        #     assert(0)
        middlebin = int(self.N_space/2)   # edges inside the source - static
        sidebin = int(middlebin/2) # edges outside the source - moving
        # dx_min = 1e-4
        # wave_source_separation = self.wave_loc_array[0,3,:] - self.x0
        # wave_sep_div = wave_source_separation / (sidebin-1)
        # index = 0
        # if wave_sep_div[index] < dx_min:
        #     index +=1
        # else:
        #     first = index
        # dx = self.wave_loc_array[0,3,index]
        dx = 1e-3
        left = np.linspace(-self.x0-dx, -self.x0, sidebin + 1)
        right = np.linspace(self.x0, self.x0 + dx, sidebin + 1)
        # middle = np.linspace(-self.x0, self.x0, middlebin + 1)
        middle = 0.5 * self.thick_quad
        self.edges = np.concatenate((left[:-1], middle[:-1], right[:])) # put them all together 
        
        # initialize derivatives
        # self.Dedges[0:sidebin] = (self.edges[0:sidebin] + self.x0 )/(self.edges[-1] - self.x0)
        # self.Dedges[sidebin:sidebin+middlebin] = 0       
        # self.Dedges[middlebin+sidebin + 1:] = (self.edges[middlebin+sidebin + 1:] - self.x0)/(self.edges[-1] - self.x0)

        self.Dedges[0] = -1.0
        self.Dedges[1] = -1.0
        self.Dedges[2] = -1.0
        self.Dedges[3:sidebin] = -(self.edges[3:sidebin] + self.x0) / (self.edges[3] + self.x0)

        # self.Dedges[int(sidebin/2+2):sidebin] = (self.edges[int(sidebin/2+2):sidebin] + self.x0 )/(self.edges[-int(sidebin/2 + 2)] - self.x0)
        self.Dedges[sidebin:sidebin+middlebin] = 0       
        self.Dedges[middlebin+sidebin + 1:] =  - np.flip(np.copy(self.Dedges[0:sidebin]))
        self.Dedges = self.Dedges * self.speed
        print(self.Dedges, 'dedges 0')
        print(self.edges, 'edges 0')


    def thick_square_init_func(self):
        print("initializing thick square source")

        dx = 1e-5

        half = int(self.N_space/2)
        self.edges = np.zeros(self.N_space+1)
        self.Dedges = np.zeros(self.N_space+1) 
        
        self.edges[half] = 0 # place center edge
        self.Dedges[half]= 0

        # self.edges[0] = -self.x0 - 2*dx# place wavefront tracking edges
        # self.edges[-1] = self.x0 + 2*dx

        # self.Dedges[0] = -1 * self.speed
        # self.Dedges[-1] = 1 * self.speed

        number_of_interior_edges = int(self.N_space/2 - 1)
        # print(number_of_interior_edges, "interior")

        # don't use N=4 
        if number_of_interior_edges == 1: # deal with N=4 case 
            self.edges[number_of_interior_edges] = -self.x0
            self.edges[number_of_interior_edges + half] = self.x0
            self.Dedges[number_of_interior_edges] = -1.0 * self.speed
            self.Dedges[number_of_interior_edges + half] = 1.0 * self.speed 
        
        else:                               # set interior edges to track the wavefront

            # set one edge to travel back towards zero and one to be fixed at the source width
            # self.set_func([half-2, half+2], [-self.x0, self.x0], [0,0])
            # self.set_func([half-1, half+1], [-self.x0+dx, self.x0-dx], [self.speed, -self.speed])
            # # set edges to track the wave
            # left_xs = np.linspace(-self.x0-2*dx, -self.x0-dx, half-2)
            # right_xs = np.linspace(self.x0+dx, self.x0+2*dx, half-2)
            # speeds = np.linspace(half-2, 1, half-2)
            # speeds = speeds/(half-2) * self.speed
            # indices_left = np.linspace(0, half-2-1, half-2)
            # indices_right = np.linspace(half+3, self.N_space, half-2)

            # self.set_func(indices_left, left_xs, -speeds)
            # self.set_func(indices_right, right_xs, np.flip(speeds))


            indices_left = np.linspace(0, half-1, half)
            indices_right = np.linspace(half+1, self.N_space, half)
            xs_left = np.zeros(half)
            xs_right = np.zeros(half)
            speeds = np.zeros(half)
            xs_left[int(half/2)] = -self.x0
            # xs_right[int(half/2)] = self.x0
            speeds[int(half/2)] = 0.0 
            xs_left[0:int(half/2)] = np.linspace(-self.x0-2*dx, -self.x0-dx,int(half/2))
            xs_left[int(half/2)+1:] = np.linspace(-self.x0+dx, -self.x0+2*dx, int(half/2)-1)
            # xs_right[0:int(half/2)] = np.linspace(self.x0-2*dx, self.x0-dx, int(half/2))
            # xs_right[int(half/2)+1:] = np.linspace(self.x0+dx, self.x0+2*dx, int(half/2)-1)
            xs_right = -np.flip(xs_left)
            speeds[0:int(half/2)] = np.linspace(int(half/2), 1, int(half/2))/int(half/2)
            speeds[int(half/2)+1:] = -np.linspace(1,int(half/2), int(half/2) -1)/ int(half/2)
            speeds = speeds * self.speed
            # print("#   #   #   #   #   #   #   #   #   #   #   ")
            # print(speeds, "speeds")
            # print("#   #   #   #   #   #   #   #   #   #   #   ")
            self.set_func(indices_left, xs_left, -speeds)
            self.set_func(indices_right, xs_right, np.flip(speeds))

            # self.edges[0:half-1] = np.linspace(-self.x0-dx, -self.x0 + dx, number_of_interior_edges + 1)
            # self.edges[half+2:] = np.linspace(self.x0 - dx, self.x0 + dx, number_of_interior_edges + 1)
            # self.edges[half-1] = -self.x0

            # # self.Dedges[1:half-1] = - self.edges[1:half-1]/self.edges[1] * self.speed
            # # self.Dedges[half+2:-1] = self.edges[half+2:-1]/self.edges[-2] * self.speed 
            # self.Dedges[0:half] = -np.linspace(1,-1, number_of_interior_edges + 1)* self.speed   
            # self.Dedges[half+1:] = np.linspace(-1,1, number_of_interior_edges + 1)* self.speed   


            self.delta_t = self.wave_loc_array[0,0,1] - self.wave_loc_array[0,0,0]
            # print(self.delta_t, 'delta_t')


    


