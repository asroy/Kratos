# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import timer_factory as timer_factory
import communicator_factory as communicator_factory
import algorithm_factory as algorithm_factory

# ==============================================================================
def CreateOptimizer( inputModelPart, optimizationSettings ):
    if  optimizationSettings["design_variables"]["variable_type"].GetString() == "vertex_morphing":
        optimizer = VertexMorphingMethod( inputModelPart, optimizationSettings )
        return optimizer
    else:
        raise NameError("Specified design control not implemented or misspelled: " + optimizationAlgorithm)

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__( self, inputModelPart, optimizationSettings ):
        
        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings
        self.__addVariablesNeededForOptimization( inputModelPart )

    # --------------------------------------------------------------------------
    def __addVariablesNeededForOptimization( self, inputModelPart ):
        inputModelPart.AddNodalSolutionStepVariable(NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        inputModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_UPDATE)
        inputModelPart.AddNodalSolutionStepVariable(DESIGN_CHANGE_ABSOLUTE)  
        inputModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        inputModelPart.AddNodalSolutionStepVariable(IS_ON_BOUNDARY)
        inputModelPart.AddNodalSolutionStepVariable(BOUNDARY_PLANE) 
        inputModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATES_DEACTIVATED) 
        inputModelPart.AddNodalSolutionStepVariable(SENSITIVITIES_DEACTIVATED) 

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.optimizationSettings["design_variables"]["input_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.inputModelPart )
        buffer_size = 1
        self.inputModelPart.SetBufferSize( buffer_size )
        self.inputModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.optimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        timer = timer_factory.CreateTimer()
        algorithmName = self.optimizationSettings["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ",timer.getTimeStamp(),": Starting optimization using the following algorithm: ", algorithmName)
        print("> ==============================================================================================================\n")
    
        designSurface = self.__getDesignSurfaceFromInputModelPart()
        listOfDampingRegions = self.__getListOfDampingRegionsFromInputModelPart()

        mapper = VertexMorphingMapper( designSurface, listOfDampingRegions, self.optimizationSettings ) 
        communicator = communicator_factory.CreateCommunicator( self.optimizationSettings )
            
        algorithm = algorithm_factory.CreateAlgorithm( designSurface, self.analyzer, mapper, communicator, self.optimizationSettings )
        algorithm.execute()       

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")
    
    # --------------------------------------------------------------------------
    def __getDesignSurfaceFromInputModelPart( self ):
        nameOfDesingSurface = self.optimizationSettings["design_variables"]["design_submodel_part_name"].GetString()
        if self.inputModelPart.HasSubModelPart( nameOfDesingSurface ):
            optimizationModel = self.inputModelPart.GetSubModelPart( nameOfDesingSurface )
            print("> The following design surface was defined:\n\n",optimizationModel)
            return optimizationModel
        else:
            raise ValueError("The following sub-model part (design surface) specified for shape optimization does not exist: ",nameOfDesingSurface)         

    # --------------------------------------------------------------------------
    def __getListOfDampingRegionsFromInputModelPart( self ):
        listOfDampingRegions = {}
        print("> The following damping regions are defined: \n")
        for regionNumber in range(self.optimizationSettings["design_variables"]["damping"]["damping_regions"].size()):
            regionName = self.optimizationSettings["design_variables"]["damping"]["damping_regions"][regionNumber]["sub_model_part_name"].GetString()
            if self.inputModelPart.HasSubModelPart(regionName):
                print(regionName)
                listOfDampingRegions[regionName] = self.inputModelPart.GetSubModelPart(regionName)
            else:
                raise ValueError("The following sub-model part specified for damping does not exist: ",regionName)         
        return listOfDampingRegions               

# ==============================================================================