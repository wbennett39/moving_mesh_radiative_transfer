import numpy as np
from numba import int64, float64
import numba
from numba.experimental import jitclass
import math
from .functions import problem_identifier

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
        # ('problem_type', int64)
        ]
#################################################################################################


    
@jitclass(data)
class mesh_class(object):
    def __init__(self, N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick):

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
        print('mesh edge velocity: ', edge_v)
        self.source_type = source_type

        if self.move_type[0] == True:
            self.move_func = 0 # simple linear
        elif self.move_type[1] == True:
            self.move_func = 1 # sqrt_t
        
        # self.problem_type = problem_identifier(self.source_typem, self.x0)
        
        self.thick = thick
        
        self.initialize_mesh()



    def move(self, t):
        """
        Called each time the rhs moves the mesh. Changes edges and Dedges
        """
        # print(self.edges)
        if self.moving == True:
            """
            This mode moves all of the edges at a constant speed,
            linearly increasing from 0 to the wavespeed
            """
            if self.move_func == 0: # simple linear
                self.edges = self.edges0 + self.Dedges*t

            elif self.move_func == 1: # sqrt_tq
                """
                This mode has the wavefront tracking edges moving at a constant speed
                and interior edges tracking the diffusive wave
        
                """
                if t > 1e-10:
                        sqrt_t = math.sqrt(t)
                else:
                    sqrt_t = math.sqrt(1e-10)
                # move the interior edges
                self.Dedges[1:-1] = self.Dedges_const[1:-1] * 2.0 * 0.5 / sqrt_t
                self.edges[1:-1] = self.edges0[1:-1] + self.Dedges_const[1:-1] * 2.0 * sqrt_t
                # move the wavefront edges
                self.Dedges[0] = self.Dedges_const[0]
                self.Dedges[-1] = self.Dedges_const[-1]
                self.edges[0] = self.edges0[0] + self.Dedges[0]*t
                self.edges[-1] = self.edges0[-1] + self.Dedges[-1]*t



    def initialize_mesh(self):
        """
        Initializes initial mesh edges and initial edge derivatives
        """
        # if self.problem_type in ['plane_IC']:
        if self.source_type[0] == 1:
            self.simple_moving_init_func()

        if self.thick == False:     # thick and thin sources have different moving functions

            # if self.problem_type in ['gaussian_IC', 'gaussian_source']:
            if self.source_type[3] == 1 or self.source_type[4] == 1:
                self.simple_moving_init_func()
            # elif self.problem_type in ['square_IC', 'square_source']:
            if self.source_type[1] == 1 or self.source_type[2] == 1:
                self.thin_square_init_func()

        elif self.thick == True:
            # if self.problem_type in ['gaussian_IC', 'gaussian_source']:
            if self.source_type[3] == 1 or self.source_type[4] == 1:
                # self.thick_gaussian_init_func()
                # have not yet optimized the mesh for this problem
                self.simple_moving_init_func()
            elif self.source_type[1] == 1 or self.source_type[2] == 1:
                self.thick_square_init_func()


        self.edges0 = self.edges
        self.Dedges_const = self.Dedges
        print(self.edges0, "edges")
        print(self.Dedges_const, "Dedges")

        if self.moving == False:
            # static mesh -- puts the edges at the final positions that the moving mesh would occupy
            # sets derivatives to 0
            self.moving = True
            self.move(self.tfinal)
            self.Dedges = self.Dedges*0
            self.moving = False

    ####### Initialization functions ########


    def simple_moving_init_func(self):
            self.edges = np.linspace(-self.x0, self.x0, self.N_space+1)
            self.Dedges = self.edges/self.edges[-1] * self.speed
    
    def thin_square_init_func(self):
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

    def thick_gaussian_init_func(self):
        return 0

    def thick_square_init_func(self):
        print("initializing thick square source")

        dx = 7

        half = int(self.N_space/2)
        self.edges = np.zeros(self.N_space+1)
        self.Dedges = np.zeros(self.N_space+1) 
        
        self.edges[half] = 0 # place center edge
        self.Dedges[half]= 0

        self.edges[0] = -self.x0 # place wavefront tracking edges
        self.edges[-1] = self.x0

        self.Dedges[0] = -self.speed
        self.Dedges[-1] = self.speed

        number_of_interior_edges = int(self.N_space/2 - 1)
        # print(number_of_interior_edges, "interior")

        if number_of_interior_edges == 1: # deal with N=4 case
            self.edges[number_of_interior_edges] = -self.x0
            self.edges[number_of_interior_edges + half] = self.x0
            self.Dedges[number_of_interior_edges] = -1.0 * self.speed
            self.Dedges[number_of_interior_edges + half] = 1.0 * self.speed 
        
        else:                               # set interior edges to track the wavefront
            self.edges[1:number_of_interior_edges+1] = np.linspace(-self.x0 + 1e-8, -self.x0 + dx/2, number_of_interior_edges )
            self.edges[half+1:-1] = np.linspace(self.x0 - dx/2, self.x0 - 1e-8, number_of_interior_edges)
            self.Dedges[1:half-1] = -1.0 * self.speed
            self.Dedges[half+2:-1] = 1.0 * self.speed         




