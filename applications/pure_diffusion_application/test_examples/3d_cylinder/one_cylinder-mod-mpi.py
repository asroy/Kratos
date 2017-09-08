##################################################################
##################################################################
#setting the domain size for the problem to be solved

domain_size = 3  # 3D problem  

#including kratos path
import sys
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.PureDiffusionApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem 
from KratosMultiphysics.ExternalSolversApplication import *
import time as time_measure

## Parse the ProjectParameters
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()
## Import KratosMPI if needed
if (parallel_type == "MPI"):
    import KratosMultiphysics.mpi as KratosMPI

### for python debugger to attach
import ptvsd
if (parallel_type == "MPI"):
    myrank = KratosMPI.mpi.rank
else:
    myrank = 0

port = 3000 + myrank

import os
print('my pid, myrank, port',os.getpid(), myrank, port,'\n')

if(myrank==0):
    ptvsd.enable_attach(secret='my_secret', address =('localhost', port))
    ptvsd.wait_for_attach()
    input()  # wait for gdb to attach


## model part definition
model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

#creating a solver object
import trilinos_pure_diffusion_solver           #we import the python file that includes the commands that we need  
solver = trilinos_pure_diffusion_solver.StaticPoissonSolverMPI(model_part,ProjectParameters)

solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
model_part.SetBufferSize(1)

# we add the DoFs  
trilinos_pure_diffusion_solver.AddDofs(model_part)


solver.Initialize()


print ("about to solve!")    
solver.Solve()
print ("Solved!")  


