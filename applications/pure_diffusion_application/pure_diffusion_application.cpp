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
#include "pure_diffusion_application.h"
#include "includes/variables.h"


namespace Kratos
{
	KRATOS_CREATE_VARIABLE(double, POINT_HEAT_SOURCE) // the other variables  needed in this app dont need to be created since they're already included in the kernel  ( CONDUCTIVITY and TEMPERATURE)

 	KratosPureDiffusionApplication::KratosPureDiffusionApplication(): //constructor  do not forget to add the ":" 
		mPoisson2D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) ),
		mPointSource2D ( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType (1) ) ) ),
		mPoissonOverset2D ( 0, Element::GeometryType::Pointer( new Line2D2 <Node<3> >( Element::GeometryType::PointsArrayType (2) ) ) ),

		mPoisson3D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),
		mPointSource3D ( 0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType (1) ) ) ),
		mPoissonOverset3D ( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType (3) ) ) )
	{}
 	
	void KratosPureDiffusionApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosPureDiffusionApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(  POINT_HEAT_SOURCE )

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Poisson2D", mPoisson2D)  //and here is our element
		KRATOS_REGISTER_CONDITION( "PointSource2D", mPointSource2D ) //and our condition
		KRATOS_REGISTER_CONDITION( "PoissonOverset2D", mPoissonOverset2D ) //and our condition

		KRATOS_REGISTER_ELEMENT("Poisson3D", mPoisson3D)  //and here is our element
		KRATOS_REGISTER_CONDITION( "PointSource3D", mPointSource3D ) //and our condition
		KRATOS_REGISTER_CONDITION( "PoissonOverset3D", mPoissonOverset3D ) //and our condition
	}


}  // namespace Kratos.