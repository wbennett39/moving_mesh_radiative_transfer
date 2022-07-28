import numpy as np
from numba import int64, float64
from numba.experimental import jitclass
import math
from .functions import problem_identifier

###############################################################################
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
        ('thick', int64)
        ]
###############################################################################

    
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
        print('mesh edge velocity: ' edge_v)
        self.source_type = source_type

        if self.move_type[0] == True:
            self.move_func = 'simple_linear'
        elif self.move_type[1] == True:
            self.move_func = 'sqrt_t'
        
        self.problem_type = problem_identifier(self.source_typem, self.x0)
        
        self.initialize_mesh()


    def move(self, t):
        """
        Called each time the rhs moves the mesh. Changes edges and Dedges
        """
        if self.moving == True:
            if self.move_func == 'simple_linear':
                self.edges = self.edges0 + self.Dedges*t
            elif self.move_func == 'sqrt_t':
                if t > 1e-10:
                        sqrt_t = math.sqrt(t)
                    else:
                        sqrt_t = math.sqrt(1e-10)
                self.Dedges = self.Dedges_const * 8.0 * 0.5 / sqrt_t
                self.edges = self.edges0 + self.Dedges_const * 8. * sqrt_t

    def initialize_mesh(self):
        """
        Initializes initial mesh edges and initial edge derivatives
        """
        if self.problem_type in ['plane_IC']:
            self.simple_moving_func()
        if self.thick == False:     # thick and thin sources have different moving functions
            if self.problem_type in ['gaussian_IC', 'gaussian_source']:
                self.simple_moving_func()
            if self.problem_type in ['square_IC', 'square_source']:
                self.thin_square_init_func()
        elif self.thick == True:
            if self.problem_type in ['gaussian_IC', 'gaussian_source']:
                self.thick_gaussian_init_func()



        self.edges0 = self.edges
        self.Dedges_const = self.Dedges

        if self.moving == False:
            # static mesh -- puts the edges at the final positions that the moving mesh would occupy
            # sets derivatives to 0
            self.move(self.tfinal)
            self.Dedges = self.Dedges*0


    def simple_moving_init_func(self):
            self.edges = np.linspace(-self.x0, self.x0, self.N_space+1)
            self.Dedges = self.edges/self.edges[-1] * self.speed
    
    def thin_square_init_func(self):

    def thick_gaussian_init_func(self):

    def thick_square_init_func(self):
