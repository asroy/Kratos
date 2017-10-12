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
#include "utilities/geometry_utilities.h" 

#include "custom_elements/PoissonElement3D.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	PoissonElement3D::PoissonElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PoissonElement3D::PoissonElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer PoissonElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new PoissonElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PoissonElement3D::~PoissonElement3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonElement3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;  //gradients matrix 
		boost::numeric::ublas::bounded_matrix<double,3,3> msD;  //conductivity matrix 
		msD = ZeroMatrix(3,3); //initializing the matrix as zero
  		array_1d<double,4> msN; //dimension = number of nodes . Position of the gauss point 
		array_1d<double,4> ms_temp; //dimension = number of nodes . . since we are using a residualbased approach 

		const unsigned int number_of_points = GetGeometry().size();

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //resizing the system in case it does not have the right size 

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Volume;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Volume); //asking for gradients and other info 

		//reading properties and conditions
		double permittivity = GetProperties()[CONDUCTIVITY];
		msD(0,0)=permittivity;
		msD(1,1)=permittivity;
		msD(2,2)=permittivity;
		
		// main loop	
		//const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

		noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));  // Bt D B

		rLeftHandSideMatrix *= Volume;


		//subtracting the dirichlet term
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp);

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonElement3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void PoissonElement3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void PoissonElement3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void PoissonElement3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}



} // Namespace Kratos
