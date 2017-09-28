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

#include "custom_conditions/PoissonOversetCondition3D.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PoissonOversetCondition3D::PoissonOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
	: 	OversetCondition(NewId, pGeometry)
{		
	//DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PoissonOversetCondition3D::PoissonOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: 	OversetCondition(NewId, pGeometry, pProperties)
{}

//************************************************************************************
//************************************************************************************
Condition::Pointer PoissonOversetCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new PoissonOversetCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

PoissonOversetCondition3D::~PoissonOversetCondition3D()
{}

//************************************************************************************
//************************************************************************************
void PoissonOversetCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	unsigned int matSize = GetGeometry().size();

	if(rRightHandSideVector.size() != matSize)
		rRightHandSideVector.resize(matSize,false);

	rRightHandSideVector = ZeroVector(matSize);

	// std::cout << __func__ <<"not implemented" << std::endl;

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void PoissonOversetCondition3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	unsigned int matSize = GetGeometry().size();

	if(rLeftHandSideMatrix.size1() != matSize)
		rLeftHandSideMatrix.resize(matSize,matSize,false);

	noalias(rLeftHandSideMatrix) = ZeroMatrix(matSize,matSize);

	if(rRightHandSideVector.size() != matSize)
		rRightHandSideVector.resize(matSize,false);

	rRightHandSideVector = ZeroVector(matSize);

	// std::cout << __func__ <<"not implemented" << std::endl;
	
	KRATOS_CATCH("")
}

void PoissonOversetCondition3D::LocalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
	std::size_t num_node = GetGeometry().size();
    rResult.resize(num_node);
    
	for (std::size_t i = 0; i < num_node; i++)
		rResult[i] = (GetGeometry()[i].GetDof(TEMPERATURE).EquationId());			
}

void PoissonOversetCondition3D::DonorEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    std::size_t num_dof = 0;

    for(std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
    {
        num_dof += rHingeDonorData(i_hinge).mEquationsId.size();
    }

    rResult.resize(num_dof);

    std::size_t i = 0;
    for(std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
    {
        const OversetAssembly::HingeDonorData & r_hinge_donor_data = rHingeDonorData(i_hinge);

        for( std::size_t j = 0; j < r_hinge_donor_data.mEquationsId.size(); j++ )
        {
            rResult[i] = r_hinge_donor_data.mEquationsId[j];
            i++;
        }
    }
}

//************************************************************************************
//************************************************************************************
void PoissonOversetCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
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
