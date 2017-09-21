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

#include "custom_conditions/PointSourceCondition3D.h"
#include "pure_diffusion_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PointSourceCondition3D::PointSourceCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PointSourceCondition3D::PointSourceCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer PointSourceCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PointSourceCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PointSourceCondition3D::~PointSourceCondition3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void PointSourceCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);

		double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
		rRightHandSideVector[0] = load;
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PointSourceCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if(rLeftHandSideMatrix.size1() != 1)
			rLeftHandSideMatrix.resize(1,1,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);
		if(rRightHandSideVector.size() != 1)
			rRightHandSideVector.resize(1,false);
		double load = GetGeometry()[0].GetSolutionStepValue(POINT_HEAT_SOURCE);
		rRightHandSideVector[0] = load;
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void PointSourceCondition3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(TEMPERATURE).EquationId());			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void PointSourceCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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
