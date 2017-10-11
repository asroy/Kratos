//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point_2d.h" 
#include "geometries/point_3d.h" 
#include "geometries/line_2d.h"
#include "includes/variables.h"

#include "pure_diffusion_application.h"

namespace Kratos
{
	KRATOS_CREATE_VARIABLE(double, POINT_HEAT_SOURCE) // the other variables  needed in this app dont need to be created since they're already included in the kernel  ( CONDUCTIVITY and TEMPERATURE)

 	KratosPureDiffusionApplication::KratosPureDiffusionApplication(): //constructor  do not forget to add the ":" 
		mPoissonElement3D          ( 0,   Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> > (   Element::GeometryType::PointsArrayType (4) ) ) ),
		mPointSourceCondition3D    ( 0, Condition::GeometryType::Pointer( new Point3D      <Node<3> > ( Condition::GeometryType::PointsArrayType (1) ) ) ),
		mPoissonHeatFluxCondition3D( 0, Condition::GeometryType::Pointer( new Point3D      <Node<3> > ( Condition::GeometryType::PointsArrayType (1) ) ) ),
		mPoissonOversetCondition3D ( 0, Condition::GeometryType::Pointer( new Triangle3D3  <Node<3> > ( Condition::GeometryType::PointsArrayType (3) ) ) )
	{}

	void KratosPureDiffusionApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosPureDiffusionApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(  POINT_HEAT_SOURCE )

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("PoissonElement3D", mPoissonElement3D)
		KRATOS_REGISTER_CONDITION( "PointSourceCondition3D", mPointSourceCondition3D )
		KRATOS_REGISTER_CONDITION( "PoissonHeatFluxCondition3D", mPoissonHeatFluxCondition3D )
		KRATOS_REGISTER_CONDITION( "PoissonOversetCondition3D", mPoissonOversetCondition3D )
	}


}  // namespace Kratos.