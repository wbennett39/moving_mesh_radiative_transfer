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
        ('debugging', int64)
        # ('problem_type', int64)
        ]
#################################################################################################


    
@jitclass(data)
class mesh_class(object):
    def __init__(self, N_space, x0, tfinal, moving, move_type, source_type, edge_v, thick):
        
        self.debugging = False

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

    def set_func(self, x_index, loc, speed):
        """
        This function takes an array of indexes and initializes the 
        edges and Dedges at the indexes given to the respective speeds and locations
        """
        for count, index in enumerate(x_index):
            self.edges[int(index)] = loc[count]
            self.Dedges[int(index)] = speed[count]


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
                self.sqrt_t_moving_func(t)
            else:
                print("no move function selected")
                assert(0)

        if self.debugging == True:
            if self.edges.all() != np.sort(self.edges).all():
                print("crossed edges")
                assert(0)

                



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
                self.thin_square_init_func()

        elif self.thick == True:
            # if self.problem_type in ['gaussian_IC', 'gaussian_source']:
            if self.source_type[3] == 1 or self.source_type[5] == 1:
                # self.thick_gaussian_init_func()
                # have not yet optimized the mesh for this problem
                self.simple_moving_init_func()

            elif self.source_type[1] == 1 or self.source_type[2] == 1:
                self.thick_square_init_func()


        self.edges0 = self.edges
        self.Dedges_const = self.Dedges


        print(self.Dedges_const, "const edges")
        print(self.edges0, 'edges0')

        # print(self.edges0, "edges")
        # print(self.Dedges_const, "Dedges")

        if self.moving == False:
            # static mesh -- puts the edges at the final positions that the moving mesh would occupy
            # sets derivatives to 0
            self.moving = True
            self.move(self.tfinal)
            self.Dedges = self.Dedges*0
            self.moving = False
            print(self.edges, "final edges")


    def sqrt_t_moving_func(self, t):
        move_factor = 3.0

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


    def thick_gaussian_init_func(self):
        return 0

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




