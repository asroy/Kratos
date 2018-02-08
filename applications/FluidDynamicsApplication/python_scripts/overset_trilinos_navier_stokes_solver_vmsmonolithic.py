from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.OversetApplication as KratosOverset       # Overset application
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD     # Fluid dynamics application

# Import base class file
import trilinos_navier_stokes_solver_vmsmonolithic

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesOversetMPISolver_VMSMonolithic(main_model_part, custom_settings)

class NavierStokesOversetMPISolver_VMSMonolithic(trilinos_navier_stokes_solver_vmsmonolithic.NavierStokesMPISolver_VMSMonolithic):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ## Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "overset_trilinos_navier_stokes_solver_vmsmonolithic",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "linear_solver_settings"       : {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-8,
                "max_levels"                         : 3,
                "symmetric"                          : false,
                "reform_preconditioner_at_each_step" : true,
                "scaling"                            : true
            },
            "volume_submodel_parts" : [""],
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "move_mesh_flag": false,
            "turbulence_model": "None"
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of NavierStokesOversetMPISolver_VMSMonolithic finished.")

    def Initialize(self):
        ## overset assembler
        self.overset_assembler = KratosOverset.OversetAssembler(self.main_model_part)

        ## add variables whose value/derivative/equation-id are needed on the hinge side
        self.overset_assembler.AddInterpolatedVariableNeedEquationId(KratosMultiphysics.VELOCITY_X)
        self.overset_assembler.AddInterpolatedVariableNeedEquationId(KratosMultiphysics.VELOCITY_Y)
        self.overset_assembler.AddInterpolatedVariableNeedEquationId(KratosMultiphysics.VELOCITY_Z)
        self.overset_assembler.AddInterpolatedVariableNeedEquationId(KratosMultiphysics.PRESSURE)

        self.overset_assembler.AddInterpolatedVariableNeedValue(KratosMultiphysics.VELOCITY_X)
        self.overset_assembler.AddInterpolatedVariableNeedValue(KratosMultiphysics.VELOCITY_Y)
        self.overset_assembler.AddInterpolatedVariableNeedValue(KratosMultiphysics.VELOCITY_Z)
        self.overset_assembler.AddInterpolatedVariableNeedValue(KratosMultiphysics.PRESSURE)

        self.overset_assembler.AddInterpolatedVariableNeedDX(KratosMultiphysics.VELOCITY_X)
        self.overset_assembler.AddInterpolatedVariableNeedDX(KratosMultiphysics.VELOCITY_Y)
        self.overset_assembler.AddInterpolatedVariableNeedDX(KratosMultiphysics.VELOCITY_Z)
        self.overset_assembler.AddInterpolatedVariableNeedDX(KratosMultiphysics.PRESSURE)

        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        ## If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        ## Creating the Trilinos convergence criteria
        self.conv_criteria = KratosTrilinos.TrilinosUPCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                               self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                               self.settings["relative_pressure_tolerance"].GetDouble(),
                                                               self.settings["absolute_pressure_tolerance"].GetDouble())

        ## Creating the Trilinos time scheme
        if (self.settings["turbulence_model"].GetString() == "None"):
            if self.settings["consider_periodic_conditions"].GetBool() == True:
                self.time_scheme = KratosOverset.OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                                self.settings["move_mesh_strategy"].GetInt(),
                                                                                                                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                                KratosCFD.PATCH_INDEX,
                                                                                                                self.overset_assembler)
            else:
                self.time_scheme = KratosOverset.OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                                self.settings["move_mesh_strategy"].GetInt(),
                                                                                                                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                                self.overset_assembler)


        ## Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        ## Construct the Trilinos builder and solver
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.builder_and_solver = KratosOverset.OversetTrilinosBlockBuilderAndSolverPeriodic(self.EpetraCommunicator,
                                                                                           guess_row_size,
                                                                                           self.trilinos_linear_solver,
                                                                                           KratosCFD.PATCH_INDEX,
                                                                                           self.overset_assembler)
        else:
            self.builder_and_solver = KratosOverset.OversetTrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                        guess_row_size,
                                                                                        self.trilinos_linear_solver,
                                                                                        self.overset_assembler)

        ## Construct the Trilinos Newton-Raphson strategy
        self.solver = KratosTrilinos.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                   self.time_scheme,
                                                                   self.trilinos_linear_solver,
                                                                   self.conv_criteria,
                                                                   self.builder_and_solver,
                                                                   self.settings["maximum_iterations"].GetInt(),
                                                                   self.settings["compute_reactions"].GetBool(),
                                                                   self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                   self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize()
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        print ("Monolithic Overset MPI solver initialization finished.")


    ## we don't want to inherit this function from the super class
    ##   because we don't want the replace_elements_and_condition_process
    def ImportModelPart(self):
        # Construct the Trilinos import model part utility
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)

        # Execute the Metis partitioning and reading
        TrilinosModelPartImporter.ExecutePartitioningAndReading()

        # Call the base class execute after reading (substitute elements, set density, viscosity and constitutie law)
        self._ExecuteAfterReading()

        # Call the base class set buffer size
        self._SetBufferSize()

        # Construct the communicators
        TrilinosModelPartImporter.CreateCommunicators()

        print ("MPI model reading finished.")


    ## we don't want to inherit this function from the super class
    ##   because we don't want the replace_elements_and_condition_process
    def _ExecuteAfterReading(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_submodel_parts",self.settings["volume_submodel_parts"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        import overset_check_and_prepare_model_process_fluid
        overset_check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

        # Read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)

        # Read the KINEMATIC VISCOSITY
        for el in self.main_model_part.Elements:
            kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            break

        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)


    def Solve(self):
        (self.solver).Solve()
