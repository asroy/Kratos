from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import pre_calculated_fluid_algorithm
BaseAlgorithm = pre_calculated_fluid_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)
        final_time = self.pp.CFD_DEM.AddEmptyValue("FinalTime").GetDouble()
        L = 0.0048 # the channel width
        center_x = 0.0044
        self.bbox_watcher = BoundingBoxRule(0.0, 2 * final_time,
                                            center_x - L, center_x + L,
                                            -0.007, -0.002,
                                            -0.005, 0.001)

    def PerformFinalOperations(self, time = None):
        BaseAlgorithm.PerformFinalOperations(self, time)
        self.particles_loader.RecordParticlesInBox(self.bbox_watcher)

    def PerformZeroStepInitializations(self):
        BaseAlgorithm.PerformZeroStepInitializations(self)
        import hdf5_io_tools
        self.particles_loader = hdf5_io_tools.ParticleHistoryLoader(self.all_model_parts.Get('SpheresPart'), self.pp, self.main_path)

    def FluidSolve(self, time = 'None', solve_system = True):
        BaseAlgorithm.FluidSolve(self, time, solve_system)
        ids = []
        X0s = []
        Y0s = []
        Z0s = []
        radii = []
        times = []

        self.disperse_phase_algorithm.watcher.GetNewParticlesData(ids, X0s, Y0s, Z0s, radii, times)
        self.particles_loader.UpdateListOfAllParticles(ids, X0s, Y0s, Z0s, radii, times)
