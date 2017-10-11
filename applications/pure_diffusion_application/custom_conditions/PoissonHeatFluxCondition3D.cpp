//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    YOUR_NAME_HERE
//

// Project includes 
#include "includes/define.h"
#include "utilities/math_utils.h"

#include "custom_conditions/PoissonHeatFluxCondition3D.h"
#include "pure_diffusion_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PoissonHeatFluxCondition3D::PoissonHeatFluxCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PoissonHeatFluxCondition3D::PoissonHeatFluxCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer PoissonHeatFluxCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PoissonHeatFluxCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PoissonHeatFluxCondition3D::~PoissonHeatFluxCondition3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonHeatFluxCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int num_node = GetGeometry().size();

		if(rRightHandSideVector.size() != num_node)
			rRightHandSideVector.resize(num_node,false);

		noalias(rRightHandSideVector) = ZeroVector(num_node);

		assert( num_node == 3 );

		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
  		array_1d<double,3> msN;
		double area;

		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, area);

		noalias(rRightHandSideVector) += ScalarVector(num_node, -1)*area;  //heat flux is -1 (heated)

		rRightHandSideVector = - rRightHandSideVector; //opposite sign please

		KRATOS_CATCH("")
	}

	void PoissonHeatFluxCondition3D::CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int num_node = GetGeometry().size();

		if( rLeftHandSideMatrix.size1() != num_node || rLeftHandSideMatrix.size2() != num_node )
			rLeftHandSideMatrix.resize(num_node,num_node,false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(num_node,num_node);

		KRATOS_CATCH("")
	}	

	//************************************************************************************
	//************************************************************************************
	void PoissonHeatFluxCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
		CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonHeatFluxCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int num_node = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(num_node*dim);
		for (int i=0;i<num_node;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(TEMPERATURE).EquationId());			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void PoissonHeatFluxCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int dim = 1;
		ConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			
			index = i*dim;
			ConditionalDofList[index] = (GetGeometry()[i].pGetDof(TEMPERATURE));
		}
	}
} // Namespace Kratos
