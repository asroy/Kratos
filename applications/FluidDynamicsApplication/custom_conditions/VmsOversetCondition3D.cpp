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

#include "custom_conditions/VmsOversetCondition3D.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
// default protected constructor
VmsOversetCondition3D::VmsOversetCondition3D(IndexType NewId)
	: 	OversetCondition(NewId)
{}

//************************************************************************************
//************************************************************************************
VmsOversetCondition3D::VmsOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
	: 	OversetCondition(NewId, pGeometry)
{}

//************************************************************************************
//************************************************************************************
VmsOversetCondition3D::VmsOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	: 	OversetCondition(NewId, pGeometry, pProperties)
{}

//************************************************************************************
//************************************************************************************
Condition::Pointer VmsOversetCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
	return Condition::Pointer(new VmsOversetCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

VmsOversetCondition3D::~VmsOversetCondition3D()
{}

//************************************************************************************
//************************************************************************************
void VmsOversetCondition3D::CalculateLocalSystem(MatrixType & rLeftHandSideMatrix, VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
{
	KRATOS_TRY

	CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
	CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void VmsOversetCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	//
	Element::GeometryType & r_parent_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	Condition::GeometryType & r_condition = GetGeometry();

	//
	const unsigned int num_dof = 4;

	unsigned int num_parent_node = r_parent_element.size();

	unsigned int size2 = num_dof*num_parent_node;
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
		size2 += num_dof*rHingeDonorData(i_hinge).NumberOfDonorNodes();

	if(rRightHandSideVector.size() != size2)
		rRightHandSideVector.resize(size2,false);

	rRightHandSideVector = ZeroVector(size2);
	
	//

	//opposite sign please
	rRightHandSideVector = - rRightHandSideVector;

	printf("rRightHandSideVector %lg %lg %lg %lg\n",rRightHandSideVector[0], rRightHandSideVector[1], rRightHandSideVector[2], rRightHandSideVector[3]);

	KRATOS_CATCH("")
}

void VmsOversetCondition3D::CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo)
{
	KRATOS_TRY

	//
	Element::GeometryType & r_parent_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	Condition::GeometryType & r_condition = GetGeometry();

	//
	const int num_dof = 4;
	unsigned int num_parent_node = r_parent_element.size();

	unsigned int size2 = num_dof*num_parent_node;
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
		size2 += num_dof*rHingeDonorData(i_hinge).NumberOfDonorNodes();

	if( rLeftHandSideMatrix.size1() != size2 || rLeftHandSideMatrix.size2() != size2 )
		rLeftHandSideMatrix.resize(size2,size2,false);

	rLeftHandSideMatrix = ZeroMatrix(size2,size2);


	KRATOS_CATCH("")
}

void VmsOversetCondition3D::LocalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

void VmsOversetCondition3D::DonorEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    std::size_t num_dof = 0;

    for(std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
    {
        num_dof += 4*rHingeDonorData(i_hinge).NumberOfDonorNodes();
    }

	rResult.clear();
    rResult.resize(num_dof);

    std::size_t i = 0;
    for(std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
    {
        const OversetAssembly::HingeDonorData & r_hinge_donor_data = rHingeDonorData(i_hinge);

		if( ! r_hinge_donor_data.IsInitialized() )
		{
			std::cout<<__func__<<"wrong! hinge donor data not Initialized"<<std::endl;
			exit(EXIT_FAILURE);
		}

        for( std::size_t j = 0; j < r_hinge_donor_data.NumberOfDonorNodes(); j++ )
        {
			rResult[i] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_X, j);
			i++;
			rResult[i] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_Y, j);
			i++;
			rResult[i] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_Z, j);
			i++;
            rResult[i] = r_hinge_donor_data.GetDonorNodeEquationId(TEMPERATURE, j);
			i++;
		}
	}

	{
		std::cout<<__func__<<": ";
		
		for( std::size_t j = 0; j < rResult.size(); j++ )
			printf(" %lu ",rResult[j]);

		std::cout<<std::endl;
	}
}

//************************************************************************************
//************************************************************************************
void VmsOversetCondition3D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
	const unsigned int num_dof = 4;

	ConditionalDofList.resize(num_dof*GetGeometry().size());

	for (unsigned int i=0; i<GetGeometry().size(); i++)
	{
		ConditionalDofList[num_dof*i]   = (GetGeometry()[i].pGetDof(VELOCITY_X));
		ConditionalDofList[num_dof*i+1] = (GetGeometry()[i].pGetDof(VELOCITY_Y));
		ConditionalDofList[num_dof*i+2] = (GetGeometry()[i].pGetDof(VELOCITY_Z));
		ConditionalDofList[num_dof*i+3] = (GetGeometry()[i].pGetDof(TEMPERATURE));
	}
}
} // Namespace Kratos
