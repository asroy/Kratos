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
#include "custom_conditions/poisson_overset_3d.h"
#include "pure_diffusion_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PoissonOverset3D::PoissonOverset3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PoissonOverset3D::PoissonOverset3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	//************************************************************************************
	//************************************************************************************
	Condition::Pointer PoissonOverset3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new PoissonOverset3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PoissonOverset3D::~PoissonOverset3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonOverset3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

    	unsigned int matSize = GetGeometry().size();

		if(rRightHandSideVector.size() != matSize)
			rRightHandSideVector.resize(matSize,false);

    	rRightHandSideVector = ZeroVector(matSize);

    	std::cout << __func__ <<"not implemented" << std::endl;

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonOverset3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

    	unsigned int matSize = GetGeometry().size();

		if(rLeftHandSideMatrix.size1() != matSize)
			rLeftHandSideMatrix.resize(matSize,matSize,false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(matSize,matSize);

		if(rRightHandSideVector.size() != matSize)
			rRightHandSideVector.resize(matSize,false);

    	rRightHandSideVector = ZeroVector(matSize);

   	 	std::cout << __func__ <<"not implemented" << std::endl;
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void PoissonOverset3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	  void PoissonOverset3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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
