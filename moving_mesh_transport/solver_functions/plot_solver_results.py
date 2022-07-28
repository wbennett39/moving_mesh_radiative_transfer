import matplotlib.pyplot as plt

from ..loading_and_saving.load_parameters import parameter_load_class

class plotter():
    def __init__(self, N_spaces, N_angles, thermal_couple, benchmarking, 
                 benchmark_xs, benchmark_phi, benchmark_e, problem_type, source):
        self.N_spaces = N_spaces
        self.N_angles = N_angles
        self.thermal_couple = thermal_couple
        self.benchmarking = benchmarking
        self.benchmark_xs = benchmark_xs
        self.benchmark_phi = benchmark_phi
        self.benchmark_e = benchmark_e
        self.problem_type = problem_type
        self.source = source
    
    def plot results(xs, phi, e_xs, e, count):
        return 0


