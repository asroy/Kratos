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
void PoissonOversetCondition3D::CalculateLocalSystem(MatrixType & rLeftHandSideMatrix, VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
{
	KRATOS_TRY

	CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
	CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

	KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void PoissonOversetCondition3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	//
	Element::GeometryType & r_parent_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	Condition::GeometryType & r_condition = GetGeometry();

	//
	unsigned int num_parent_node = r_parent_element.size();

	unsigned int size2 = num_parent_node;
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
		size2 += rHingeDonorData(i_hinge).mNs.size();

	if(rRightHandSideVector.size() != size2)
		rRightHandSideVector.resize(size2,false);

	rRightHandSideVector = ZeroVector(size2);

	//
	Vector parent_element_nodes_temp(r_parent_element.size());
	for(std::size_t i = 0; i < r_parent_element.size(); i++)
		parent_element_nodes_temp[i] = r_parent_element[i].FastGetSolutionStepValue(TEMPERATURE);

	printf("%s %d parent_element_nodes_temp %lg %lg %lg %lg\n",__func__,1, parent_element_nodes_temp[0], parent_element_nodes_temp[1], parent_element_nodes_temp[2], parent_element_nodes_temp[3]);
		
	//loop over hinges
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
	{
		const IntegrationPointType & r_condition_integration_point = rHinge(i_hinge);
		const double weight = r_condition_integration_point.Weight();
		const OversetAssembly::HingeDonorData & r_hinge_donor_data = rHingeDonorData(i_hinge);

		//condition
		//  condition shape functions
		Vector Ns_condition;
		r_condition.ShapeFunctionsValues( Ns_condition, r_condition_integration_point );

		//  condition jacobian
		Matrix condition_jacobian;
		r_condition.Jacobian(condition_jacobian, r_condition_integration_point);
		const double condition_jacobian_determinant = r_condition.DeterminantOfJacobian(r_condition_integration_point);

		//  condition normal vector
		Vector normal_vector = ConditionJacobianToOutwardNormalVector(condition_jacobian);

		//parent element
		//  parent element point
		Point parent_element_point;
		r_parent_element.PointLocalCoordinates( parent_element_point, HingeGlobalCoordinate(i_hinge) );

		//  parent element shape functions
		Vector Ns_parent_element(r_parent_element.size());
		Matrix DNs_DEs_parent_element( r_parent_element.size(), 3 );
		r_parent_element.ShapeFunctionsValues( Ns_parent_element, parent_element_point );
		r_parent_element.ShapeFunctionsLocalGradients(DNs_DEs_parent_element, parent_element_point);

		printf("%s %d Ns_parent_element %lg %lg %lg %lg\n",__func__,21, Ns_parent_element[0], Ns_parent_element[1], Ns_parent_element[2], Ns_parent_element[3]);
		

		//  parent element Jacbian matrix
		Matrix jinv_parent_element(3,3);
		r_parent_element.InverseOfJacobian( jinv_parent_element, parent_element_point );

		//  parent shape functions Gradient
		Matrix DNs_DXs_parent_element = prod(DNs_DEs_parent_element, jinv_parent_element);

		//  parent element temperature value and global gradient
		double temp_parent_element = boost::numeric::ublas::inner_prod( Ns_parent_element, parent_element_nodes_temp );
		
		Vector Dtemp_DXs_parent_element = prod( trans(DNs_DXs_parent_element), parent_element_nodes_temp );

		printf("%s %d temp_parent_element %lg Dtemp_DXs_parent_element %lg %lg %lg\n",__func__,22, temp_parent_element, Dtemp_DXs_parent_element[0], Dtemp_DXs_parent_element[1], Dtemp_DXs_parent_element[2]);
		
		//donor element
		//  donor element temperature and global gradient
		double temp_donor_element = r_hinge_donor_data.GetValue(TEMPERATURE);
		std::vector<double> dtemp_dxs = r_hinge_donor_data.GetDVDXs(TEMPERATURE);

		Vector Dtemp_DXs_donor_element(3);
		Dtemp_DXs_donor_element[0] = dtemp_dxs[0];
		Dtemp_DXs_donor_element[1] = dtemp_dxs[1];
		Dtemp_DXs_donor_element[2] = dtemp_dxs[2];

		printf("%s %d temp_donor_element %lg Dtemp_DXs_donor_element %lg %lg %lg\n",__func__,23, temp_donor_element, Dtemp_DXs_donor_element[0], Dtemp_DXs_donor_element[1], Dtemp_DXs_donor_element[2]);

		//calculate mu
		double mu = 0;
		for( std::size_t i = 0; i < num_parent_node; i++ )
			mu += std::abs(DNs_DXs_parent_element(i,0)) + std::abs(DNs_DXs_parent_element(i,1)) + std::abs(DNs_DXs_parent_element(i,2));

		mu /= 3*num_parent_node;

		mu *= 1000;

		printf("%s %d weight %lg determ %lg normal_vector %lg %lg %lg\n",__func__,24, weight, condition_jacobian_determinant, normal_vector[0], normal_vector[1], normal_vector[2]);

		//calculate rhs
		for( std::size_t i_parent_element_node = 0; i_parent_element_node < r_parent_element.size(); i_parent_element_node++ )
		{
			MatrixRow DN_DXs_parent_element_i(DNs_DXs_parent_element,i_parent_element_node);

			rRightHandSideVector[i_parent_element_node] += weight * condition_jacobian_determinant *
				( 
					- 0.5 * boost::numeric::ublas::inner_prod( ( Dtemp_DXs_parent_element + Dtemp_DXs_donor_element ), normal_vector ) * Ns_parent_element[i_parent_element_node]
					+ 0.5 * ( temp_parent_element - temp_donor_element ) * boost::numeric::ublas::inner_prod( normal_vector, DN_DXs_parent_element_i )
					+ mu * ( temp_parent_element - temp_donor_element ) * Ns_parent_element[i_parent_element_node]
				);
		}
	}

	//opposite sign please
	rRightHandSideVector = - rRightHandSideVector;

	printf("rRightHandSideVector %lg %lg %lg %lg\n",rRightHandSideVector[0], rRightHandSideVector[1], rRightHandSideVector[2], rRightHandSideVector[3]);

	KRATOS_CATCH("")
}

void PoissonOversetCondition3D::CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo)
{
	KRATOS_TRY

	//
	Element::GeometryType & r_parent_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	Condition::GeometryType & r_condition = GetGeometry();

	//
	unsigned int num_parent_node = r_parent_element.size();

	unsigned int size2 = num_parent_node;
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
		size2 += rHingeDonorData(i_hinge).mNs.size();

	if( rLeftHandSideMatrix.size1() != size2 || rLeftHandSideMatrix.size2() != size2 )
		rLeftHandSideMatrix.resize(size2,size2,false);

	rLeftHandSideMatrix = ZeroMatrix(size2,size2);

	//
	Vector parent_element_nodes_temp(r_parent_element.size());
	for(std::size_t i = 0; i < r_parent_element.size(); i++)
		parent_element_nodes_temp[i] = r_parent_element[i].FastGetSolutionStepValue(TEMPERATURE);

	//loop over hinges
	std::size_t j_dof_base = num_parent_node;

	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
	{
		const IntegrationPointType & r_condition_integration_point = rHinge(i_hinge);
		const double weight = r_condition_integration_point.Weight();
		const OversetAssembly::HingeDonorData & r_hinge_donor_data = rHingeDonorData(i_hinge);

		//condition
		//  condition shape functions
		Vector Ns_condition;
		r_condition.ShapeFunctionsValues( Ns_condition, r_condition_integration_point );

		//  condition jacobian
		Matrix condition_jacobian;
		r_condition.Jacobian(condition_jacobian, r_condition_integration_point);
		const double condition_jacobian_determinant = r_condition.DeterminantOfJacobian(r_condition_integration_point);

		//  condition normal vector
		Vector normal_vector = ConditionJacobianToOutwardNormalVector(condition_jacobian);

		//parent element
		//  parent element point
		Point parent_element_point;
		r_parent_element.PointLocalCoordinates( parent_element_point, HingeGlobalCoordinate(i_hinge) );

		//  parent element shape functions
		Vector Ns_parent_element(r_parent_element.size());
		Matrix DNs_DEs_parent_element( r_parent_element.size(), 3 );
		r_parent_element.ShapeFunctionsValues( Ns_parent_element, parent_element_point );
		r_parent_element.ShapeFunctionsLocalGradients(DNs_DEs_parent_element, parent_element_point);

		//  parent element Jacbian matrix
		Matrix jinv_parent_element(3,3);
		r_parent_element.InverseOfJacobian( jinv_parent_element, parent_element_point );

		//  parent shape functions Gradient
		Matrix DNs_DXs_parent_element = prod(DNs_DEs_parent_element, jinv_parent_element);

		//donor element
		std::size_t num_donor_node = r_hinge_donor_data.NumberOfDonorNodes();

		//  donor element shape function value adn global gradient
		Vector Ns_donor_element(num_donor_node);
		for( std::size_t i = 0; i < num_donor_node; i++ )
			Ns_donor_element[i] = r_hinge_donor_data.mNs[i];

		Matrix DNs_DXs_donor_element(num_donor_node,3);
		for( std::size_t i = 0; i < num_donor_node; i++ )
		{
			DNs_DXs_donor_element(i,0) = r_hinge_donor_data.mDNsDXs[i][0];
			DNs_DXs_donor_element(i,1) = r_hinge_donor_data.mDNsDXs[i][1];
			DNs_DXs_donor_element(i,2) = r_hinge_donor_data.mDNsDXs[i][2];
		}

		//calculate mu
		double mu = 0;
		for( std::size_t i = 0; i < num_parent_node; i++ )
			mu += std::abs(DNs_DXs_parent_element(i,0)) + std::abs(DNs_DXs_parent_element(i,1)) + std::abs(DNs_DXs_parent_element(i,2));
		
		mu /= 3*num_parent_node;

		mu *= 1000;

		//calculate lhs
		for( std::size_t i_parent_element_node = 0; i_parent_element_node < num_parent_node; i_parent_element_node++ )
		{
			MatrixRow DN_DXs_parent_element_i(DNs_DXs_parent_element,i_parent_element_node);

			for( std::size_t j_parent_element_node = 0; j_parent_element_node < num_parent_node; j_parent_element_node++ )
			{
				MatrixRow DN_DXs_parent_element_j(DNs_DXs_parent_element,j_parent_element_node);
	
				rLeftHandSideMatrix(i_parent_element_node,j_parent_element_node) += weight * condition_jacobian_determinant *
					(
						- 0.5 * boost::numeric::ublas::inner_prod( DN_DXs_parent_element_j, normal_vector ) * Ns_parent_element[i_parent_element_node]
						+ 0.5 * Ns_parent_element[j_parent_element_node] * boost::numeric::ublas::inner_prod( normal_vector, DN_DXs_parent_element_i )
						+ mu * Ns_parent_element[j_parent_element_node] * Ns_parent_element[i_parent_element_node]
					);
			}

			for( std::size_t j_donor_element_node = 0; j_donor_element_node < num_donor_node; j_donor_element_node++ )
			{
				MatrixRow DN_DXs_donor_element_j(DNs_DXs_donor_element,j_donor_element_node);
	
				rLeftHandSideMatrix(i_parent_element_node,j_dof_base + j_donor_element_node) += weight * condition_jacobian_determinant *
					( 
						- 0.5 * boost::numeric::ublas::inner_prod( DN_DXs_donor_element_j, normal_vector ) * Ns_parent_element[i_parent_element_node]
						- 0.5 * Ns_donor_element[j_donor_element_node] * boost::numeric::ublas::inner_prod( normal_vector, DN_DXs_parent_element_i )
						- mu * Ns_donor_element[j_donor_element_node] * Ns_parent_element[i_parent_element_node]
					);
			}
		}

		j_dof_base += num_donor_node;
	}

	KRATOS_CATCH("")
}

void PoissonOversetCondition3D::LocalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
	Element::GeometryType & r_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	
	std::size_t num_node = r_element.size();

	rResult.clear();
    rResult.resize(num_node);
    
	for (std::size_t i = 0; i < num_node; i++)
		rResult[i] = r_element[i].GetDof(TEMPERATURE).EquationId();

	{
		std::cout<<__func__<<": ";
		
		for( std::size_t j = 0; j < rResult.size(); j++ )
			printf(" %lu ",rResult[j]);

		std::cout<<std::endl;
	}
}

void PoissonOversetCondition3D::DonorEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    std::size_t num_dof = 0;

    for(std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
    {
        num_dof += rHingeDonorData(i_hinge).NumberOfDonorNodes();
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
