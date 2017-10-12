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

		const double const_heat_flux = -1; //heat goes in

		Condition::GeometryType & r_condition = GetGeometry();
		
		const unsigned int num_condition_node = r_condition.size();

		if(rRightHandSideVector.size() != num_condition_node)
			rRightHandSideVector.resize(num_condition_node,false);

		noalias(rRightHandSideVector) = ZeroVector(num_condition_node);

		const IntegrationPointsArrayType & r_condition_integration_points = r_condition.IntegrationPoints();
		
		//loop over hinges
		for( const IntegrationPointType & r_condition_integration_point : r_condition_integration_points )
		{
			const double weight = r_condition_integration_point.Weight();

			//condition
			//  condition shape functions
			Vector Ns_condition;
			r_condition.ShapeFunctionsValues( Ns_condition, r_condition_integration_point );

			//  condition jacobian
			const double condition_jacobian_determinant = r_condition.DeterminantOfJacobian(r_condition_integration_point);

			//calculate rhs
			for( std::size_t i_condition_node = 0; i_condition_node < num_condition_node; i_condition_node++ )
				rRightHandSideVector[i_condition_node] += Ns_condition[i_condition_node] * weight * condition_jacobian_determinant * const_heat_flux;
		}

		rRightHandSideVector = - rRightHandSideVector; //opposite sign please

		KRATOS_CATCH("")
	}

	void PoissonHeatFluxCondition3D::CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int num_condition_node = GetGeometry().size();

		if( rLeftHandSideMatrix.size1() != num_condition_node || rLeftHandSideMatrix.size2() != num_condition_node )
			rLeftHandSideMatrix.resize(num_condition_node,num_condition_node,false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(num_condition_node,num_condition_node);

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
		int num_condition_node = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 1;
		rResult.resize(num_condition_node*dim);
		for (int i=0;i<num_condition_node;i++)
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
