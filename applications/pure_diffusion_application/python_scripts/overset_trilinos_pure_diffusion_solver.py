#importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.OversetApplication as KratosOverset
import KratosMultiphysics.PureDiffusionApplication as KratosDiffuson

def AddDofs(model_part):
    for node in model_part.Nodes:

        #adding dofs
        node.AddDof(KratosMultiphysics.TEMPERATURE);

    print ("variables for the Poisson solver added correctly")

class StaticPoissonSolverMPI:
    #######################################################################
    def __init__(self,model_part,custom_settings):  #constructor of the class 
        self.model_part = model_part

        default_settings = KratosMultiphysics.Parameters("""
        {
            "problem_data"                     : {
                "problem_name"    : "unknown",
                "model_part_name" : "unknown",
                "domain_size"     : 3,
                "parallel_type"   : "unknown"
            },
            "solver_settings"                  : {
                "model_import_settings"        : {
                    "input_type"     : "mpda",
                    "input_filename" : "unknown"
                }
            },
            "linear_solver_settings"       : {
                "solver_type"                        : "MultiLevelSolver",
                "max_iteration"                      : 200,
                "tolerance"                          : 1e-8,
                "max_levels"                         : 3,
                "symmetric"                          : false,
                "reform_preconditioner_at_each_step" : true,
                "scaling"                            : true
            }
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
    #######################################################################
    def Initialize(self):

        #overset assembler
        self.overset_assembler = KratosOverset.OversetAssembler(self.model_part)

        #definition of the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])  #we set the type of solver that we want 
        
        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Creating the Trilinos time scheme
        self.time_scheme = KratosOverset.OversetTrilinosResidualBasedIncrementalUpdateStaticScheme(self.overset_assembler)

        ## Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver
        if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        ## Construct the Trilinos builder and solver
        self.builder_and_solver = KratosOverset.OversetTrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                    guess_row_size,
                                                                                    self.trilinos_linear_solver,
                                                                                    self.overset_assembler)


        ## Construct strategy
        #option 1: using linear strategy
        # CalculateReactionFlag = False
        # ReformDofSetAtEachStep = False
        # CalculateNormDxFlag = False
        # MoveMeshFlag = False
        # self.solver = KratosTrilinos.TrilinosLinearStrategy(self.model_part,
        #                                                     self.time_scheme,
        #                                                     self.trilinos_linear_solver,
        #                                                     self.builder_and_solver,
        #                                                     CalculateReactionFlag,
        #                                                     ReformDofSetAtEachStep,
        #                                                     CalculateNormDxFlag,
        #                                                     MoveMeshFlag)

        #option 2: using newton method
        self.convergence_criterion = KratosTrilinos.TrilinosResidualCriteria(1e-10, 1e-14)
        MaxNewtonIteration = 30
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        self.solver = KratosTrilinos.TrilinosNewtonRaphsonStrategy(self.model_part,
                                                            self.time_scheme,
                                                            self.trilinos_linear_solver,
                                                            self.convergence_criterion,
                                                            self.builder_and_solver,
                                                            MaxNewtonIteration,
                                                            CalculateReactionFlag,
                                                            ReformDofSetAtEachStep,
                                                            MoveMeshFlag)

        echo_level = 0
        (self.solver).SetEchoLevel(echo_level)

        # self.solver.Initialize()
        self.solver.Check()

        print ("MPI solver initialization finished.")

    #######################################################################
    def AddVariables(self):
        ## Add variables from the base class
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE);
        self.model_part.AddNodalSolutionStepVariable(KratosDiffuson.POINT_HEAT_SOURCE);

        ## Add specific MPI variables
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMPI.mpi.world.barrier()

        if KratosMPI.mpi.rank == 0:
            print("Variables added correctly in each processor.")

    #######################################################################
    def ImportModelPart(self):
        # Construct the Trilinos import model part utility
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.model_part, self.settings["solver_settings"])

        # Execute the Metis partitioning and reading
        TrilinosModelPartImporter.ExecutePartitioningAndReading()

        # Construct the communicators
        TrilinosModelPartImporter.CreateCommunicators()

        print ("MPI model reading finished.")
                 
    #######################################################################   
    def Solve(self):
        self.overset_assembler.GenerateHinges()
        self.overset_assembler.SearchHingesDonor()

        print("overset_trilinos_pure_diffusion_solver::Solve")
        
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
