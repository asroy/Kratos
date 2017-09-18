//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "overset_application.h"
#include "includes/variables.h"


namespace Kratos
{
   	KRATOS_CREATE_VARIABLE(int, BLOCK_ID)

 	KratosOversetApplication::KratosOversetApplication()
 	{}
 	
 	void KratosOversetApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosOversetApplication... " << std::endl;
 
		KRATOS_REGISTER_VARIABLE( BLOCK_ID )

		KRATOS_REGISTER_CONDITION( "OversetCondition", mOversetCondition )

 	}

}  // namespace Kratos.