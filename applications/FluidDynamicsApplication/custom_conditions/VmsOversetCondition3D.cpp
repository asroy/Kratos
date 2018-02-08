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
	VectorType parent_element_nodes_velocity_x(r_parent_element.size());
	VectorType parent_element_nodes_velocity_y(r_parent_element.size());
	VectorType parent_element_nodes_velocity_z(r_parent_element.size());
	VectorType parent_element_nodes_pressure  (r_parent_element.size());

	for(std::size_t i = 0; i < r_parent_element.size(); i++)
	{
		parent_element_nodes_velocity_x[i] = r_parent_element[i].FastGetSolutionStepValue(VELOCITY_X);
		parent_element_nodes_velocity_y[i] = r_parent_element[i].FastGetSolutionStepValue(VELOCITY_Y);
		parent_element_nodes_velocity_z[i] = r_parent_element[i].FastGetSolutionStepValue(VELOCITY_Z);
		parent_element_nodes_pressure  [i] = r_parent_element[i].FastGetSolutionStepValue(PRESSURE);
	}

	printf("%s %d parent_element_nodes_pressure %lg %lg %lg %lg\n",__func__,1, parent_element_nodes_pressure[0], parent_element_nodes_pressure[1], parent_element_nodes_pressure[2], parent_element_nodes_pressure[3]);

	//loop over hinges
	for( std::size_t i_hinge = 0; i_hinge < NumberOfHinges(); i_hinge++ )
	{
		const IntegrationPointType & r_condition_integration_point = rHinge(i_hinge);
		const double weight = r_condition_integration_point.Weight();
		const OversetAssembly::HingeDonorData & r_hinge_donor_data = rHingeDonorData(i_hinge);

		//condition
		//  condition shape functions
		VectorType Ns_condition;
		r_condition.ShapeFunctionsValues( Ns_condition, r_condition_integration_point );

		//  condition jacobian
		Matrix condition_jacobian;
		r_condition.Jacobian(condition_jacobian, r_condition_integration_point);
		const double condition_jacobian_determinant = r_condition.DeterminantOfJacobian(r_condition_integration_point);

		//  condition normal vector
		VectorType normal_vector = ConditionJacobianToOutwardNormalVector(condition_jacobian);

		//parent element
		//  parent element point
		Point parent_element_point;
		r_parent_element.PointLocalCoordinates( parent_element_point, HingeGlobalCoordinate(i_hinge) );

		//  parent element shape functions
		VectorType Ns_parent_element(r_parent_element.size());
		Matrix DNs_DEs_parent_element( r_parent_element.size(), 3 );
		r_parent_element.ShapeFunctionsValues( Ns_parent_element, parent_element_point );
		r_parent_element.ShapeFunctionsLocalGradients(DNs_DEs_parent_element, parent_element_point);

		printf("%s %d Ns_parent_element %lg %lg %lg %lg\n",__func__,21, Ns_parent_element[0], Ns_parent_element[1], Ns_parent_element[2], Ns_parent_element[3]);
		
		//  parent element Jacbian matrix
		Matrix jinv_parent_element(3,3);
		r_parent_element.InverseOfJacobian( jinv_parent_element, parent_element_point );

		//  parent shape functions Gradient
		Matrix DNs_DXs_parent_element = prod(DNs_DEs_parent_element, jinv_parent_element);

		//  parent element value and global gradient
		VectorType parent_solution(4);
		parent_solution[0] = boost::numeric::ublas::inner_prod( Ns_parent_element, parent_element_nodes_velocity_x );
		parent_solution[1] = boost::numeric::ublas::inner_prod( Ns_parent_element, parent_element_nodes_velocity_y );
		parent_solution[2] = boost::numeric::ublas::inner_prod( Ns_parent_element, parent_element_nodes_velocity_z );
		parent_solution[3] = boost::numeric::ublas::inner_prod( Ns_parent_element, parent_element_nodes_pressure );
		
		printf("%s %d parent_element_pressure %lg\n",__func__,22, parent_solution[3]);

		//   parent flux operator
		MatrixType parent_flux_operator;
		CalculateFluxOperator( parent_flux_operator, parent_solution, normal_vector );

		VectorType parent_solution_flux = prod( parent_flux_operator, parent_solution );

		//donor element
		//  donor element temperature and global gradient
		VectorType donor_solution(4);
		donor_solution[0] = r_hinge_donor_data.GetValue(VELOCITY_X);
		donor_solution[1] = r_hinge_donor_data.GetValue(VELOCITY_Y);
		donor_solution[2] = r_hinge_donor_data.GetValue(VELOCITY_Z);
		donor_solution[3] = r_hinge_donor_data.GetValue(PRESSURE);

		printf("%s %d donor_element_pressure %lg \n",__func__,23, donor_solution[3]);

		//   donor flux operator
		MatrixType donor_flux_operator;
		CalculateFluxOperator( donor_flux_operator, donor_solution, normal_vector );

		VectorType donor_solution_flux = prod( donor_flux_operator, donor_solution );

		//calculate mu
		MatrixType mu_matrix;
		CalculateMuMatrix( mu_matrix, parent_solution, donor_solution, normal_vector );

		mu_matrix *= 1000;

		printf("%s %d weight %lg determ %lg normal_vector %lg %lg %lg\n",__func__,24, weight, condition_jacobian_determinant, normal_vector[0], normal_vector[1], normal_vector[2]);

		//calculate rhs
		unsigned int i_dof = 0;
		for( std::size_t i_parent_element_node = 0; i_parent_element_node < r_parent_element.size(); i_parent_element_node++ )
		{
			double tmp = Ns_parent_element[i_parent_element_node];
			VectorType parent_test_function(4);
			parent_test_function[0] = tmp;
			parent_test_function[1] = tmp;
			parent_test_function[2] = tmp;
			parent_test_function[3] = tmp;

			VectorType parent_test_function_flux = prod( parent_flux_operator, parent_test_function );

			VectorType mm = prod( mu_matrix, parent_test_function );
	
			for( int ii = 0; ii < 4; ii++ )
			{
				rRightHandSideVector[i_dof++] += weight * condition_jacobian_determinant *
					( 
						- 0.5 * ( parent_solution_flux[ii] + donor_solution_flux[ii] ) * parent_test_function[ii]
						+ 0.5 * ( parent_solution[ii] - donor_solution[ii] ) * parent_test_function_flux[ii]
						+ ( parent_solution[ii] - donor_solution[ii] ) * mm[ii]
					);
			}
		}
	}

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
	Element::GeometryType & r_element = ( const_cast<Element &> (rAdjacentElement()) ).GetGeometry();
	
	std::size_t num_node = r_element.size();

    const SizeType num_dof = 4*num_node;
    unsigned int i = 0;

    if (rResult.size() != num_dof)
        rResult.resize(num_dof, false);

    for (unsigned int iNode = 0; iNode < num_node; ++iNode)
    {
        rResult[i++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[i++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[i++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[i++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
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
			rResult[i++] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_X, j);
			rResult[i++] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_Y, j);
			rResult[i++] = r_hinge_donor_data.GetDonorNodeEquationId(VELOCITY_Z, j);
            rResult[i++] = r_hinge_donor_data.GetDonorNodeEquationId(PRESSURE, j);
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

void VmsOversetCondition3D::CalculateFluxOperator( MatrixType & r_flux_operator, const VectorType & r_solution, const VectorType & r_normal )
{
	r_flux_operator = ZeroMatrix(4,4);

	double div = r_solution[0]*r_normal[0] + r_solution[1]*r_normal[1] + r_solution[2]*r_normal[2];

	r_flux_operator(0,0) = div;
	r_flux_operator(1,1) = div;
	r_flux_operator(2,2) = div;

	r_flux_operator(0,3) = - r_normal(0);
	r_flux_operator(1,3) = - r_normal(1);
	r_flux_operator(2,3) = - r_normal(2);

	r_flux_operator(3,0) = r_normal(0);
	r_flux_operator(3,1) = r_normal(1);
	r_flux_operator(3,2) = r_normal(2);
}

void VmsOversetCondition3D::CalculateMuMatrix( MatrixType & r_mu_matrix, const VectorType & r_parent_solution, const VectorType & r_donor_solution, const VectorType & r_normal )
{
	r_mu_matrix = ZeroMatrix(4,4);

	double parent_div = std::abs(r_parent_solution[0]*r_normal[0]) + std::abs(r_parent_solution[1]*r_normal[1]) + std::abs(r_parent_solution[2]*r_normal[2]);
	double donor_div = std::abs(r_donor_solution[0]*r_normal[0]) + std::abs(r_donor_solution[1]*r_normal[1]) + std::abs(r_donor_solution[2]*r_normal[2]);

	double div = 0.5*( parent_div + donor_div );
	double n0 = std::abs(r_normal[0]);
	double n1 = std::abs(r_normal[1]);
	double n2 = std::abs(r_normal[2]);

	r_mu_matrix(0,0) = div;
	r_mu_matrix(1,1) = div;
	r_mu_matrix(2,2) = div;

	r_mu_matrix(0,3) = -n0;
	r_mu_matrix(1,3) = -n1;
	r_mu_matrix(2,3) = -n2;

	r_mu_matrix(3,0) = n0;
	r_mu_matrix(3,1) = n1;
	r_mu_matrix(3,2) = n2;
}

} // Namespace Kratos
