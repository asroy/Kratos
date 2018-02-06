/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 


// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "linear_solvers/linear_solver.h"

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "trilinos_application/trilinos_application.h"
#include "trilinos_application/trilinos_space.h"

//overset includes
#include "custom_python/add_overset_strategies_to_python.h"
#include "custom_strategies/builder_and_solvers/overset_trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/overset_trilinos_block_builder_and_solver_periodic.h"

namespace Kratos
{
	namespace Python
	{		
		using namespace boost::python;

		void  AddOversetStrategiesToPython()
		{
			typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
			typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
		
			using TrilinosBlockBuilderAndSolverType = TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
	
			using OversetTrilinosBlockBuilderAndSolverType = OversetTrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;

			using TrilinosBlockBuilderAndSolverPeriodicType = TrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
	
			using OversetTrilinosBlockBuilderAndSolverPeriodicType = OversetTrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;

			class_< OversetTrilinosBlockBuilderAndSolverType,bases< TrilinosBlockBuilderAndSolverType >, boost::noncopyable >
			( "OversetTrilinosBlockBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, OversetAssembly::OversetAssembler & > () );

			class_< OversetTrilinosBlockBuilderAndSolverPeriodicType,bases< TrilinosBlockBuilderAndSolverPeriodicType >, boost::noncopyable >
			( "OversetTrilinosBlockBuilderAndSolverPeriodic", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, OversetAssembly::OversetAssembler & > () );
		}

	}  // namespace Python.

} // Namespace Kratos

