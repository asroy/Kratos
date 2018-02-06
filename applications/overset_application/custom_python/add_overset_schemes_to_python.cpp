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
#include "custom_python/add_overset_schemes_to_python.h"
#include "custom_strategies/schemes/overset_trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/schemes/overset_trilinos_predictorcorrector_velocity_bossak_scheme_turbulent.h"

namespace Kratos
{
	namespace Python
	{		
		using namespace boost::python;

		void  AddOversetSchemesToPython()
		{
			typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
			typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

			using TrilinosResidualBasedIncrementalUpdateStaticSchemeType = TrilinosResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType> ;
		
			using OversetTrilinosResidualBasedIncrementalUpdateStaticSchemeType = OversetTrilinosResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType> ;

			class_< OversetTrilinosResidualBasedIncrementalUpdateStaticSchemeType,bases< TrilinosResidualBasedIncrementalUpdateStaticSchemeType >, boost::noncopyable >
			( "OversetTrilinosResidualBasedIncrementalUpdateStaticScheme", init<OversetAssembly::OversetAssembler & > () );


			using TrilinosPredictorCorrectorVelocityBossakSchemeTurbulentType = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType> ;

			using OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulentType = OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType> ;
			
			class_< OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulentType,bases< TrilinosPredictorCorrectorVelocityBossakSchemeTurbulentType >, boost::noncopyable >
			( "OversetTrilinosPredictorCorrectorVelocityBossakSchemeTurbulent", init<OversetAssembly::OversetAssembler & > () );
		}

	}  // namespace Python.

} // Namespace Kratos

