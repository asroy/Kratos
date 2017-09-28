//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    Chao Liu
//

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"

#include "overset_application.h"
#include "custom_conditions/OversetCondition.h"

namespace Kratos
{
//default protected constructor
OversetCondition::OversetCondition(IndexType NewId)
    :   Condition(NewId),
        mpAdjacentElement{nullptr},
        mAdjacentElementSide{0},
        mIntegrationMethod{GetIntegrationMethod()}
{}

OversetCondition::OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    :   Condition(NewId, pGeometry),
        mpAdjacentElement{nullptr},
        mAdjacentElementSide{0},
        mIntegrationMethod{GetIntegrationMethod()}
{}

OversetCondition::OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :   Condition(NewId, pGeometry, pProperties),
        mpAdjacentElement{nullptr},
        mAdjacentElementSide{0},
        mIntegrationMethod{GetIntegrationMethod()}
{}

OversetCondition::~OversetCondition()
{}

Condition::Pointer OversetCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new OversetCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void OversetCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    EquationIdVectorType local_equations_id;
    EquationIdVectorType donor_equations_id;

    LocalEquationIdVector(local_equations_id, rCurrentProcessInfo);
    DonorEquationIdVector(donor_equations_id, rCurrentProcessInfo);

    rResult.resize(local_equations_id.size()+donor_equations_id.size());
    
    std::size_t i = 0;

    for( std::size_t j = 0; j < local_equations_id.size(); j++ )
    {
        rResult[i] = local_equations_id[j];
        i++;
    }

    for( std::size_t j = 0; j < donor_equations_id.size(); j++ )
    {
        rResult[i] = donor_equations_id[j];
        i++;
    }
}

void OversetCondition::LocalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    KRATOS_ERROR << "Calling base class 'LocalEquationIdVector' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
}

void OversetCondition::DonorEquationIdVector(EquationIdVectorType& rResult, ProcessInfo & rCurrentProcessInfo)
{
    KRATOS_ERROR << "Calling base class 'DonorEquationIdVector' method instead of derived class one. Please check the definition of derived class. " << *this << std::endl;
}

void OversetCondition::GenerateHinges()
{
    mHinges.clear();
    mHingesDonorInfo.clear();
    mHingesDonorData.clear();

    const IntegrationPointsArrayType & r_integration_points = GetGeometry().IntegrationPoints(mIntegrationMethod);
    const std::size_t num_hinge = r_integration_points.size();

    //mHinges
    mHinges.reserve(num_hinge);
    for( const OversetCondition::IntegrationPointType r_integration_point : r_integration_points )
        mHinges.push_back(r_integration_point);

    //mHingesDonorInfo, mHingesDonorData
    mHingesDonorInfo.resize(num_hinge);
    mHingesDonorData.resize(num_hinge);
}

const std::size_t OversetCondition::NumberOfHinges() const
{ return mHinges.size(); }

const OversetCondition::IntegrationPointType & OversetCondition::rHinge(const std::size_t i_hinge) const
{ return mHinges[i_hinge]; }

OversetAssembly::HingeDonorInfo & OversetCondition::rHingeDonorInfo(const std::size_t i_hinge)
{ return mHingesDonorInfo[i_hinge]; }

OversetAssembly::HingeDonorData & OversetCondition::rHingeDonorData(const std::size_t i_hinge)
{ return mHingesDonorData[i_hinge]; }

const OversetCondition::PointType OversetCondition::HingeGlobalCoordinate(const std::size_t i_hinge ) const
{
    const GeometryType & r_nodes = GetGeometry();

    Vector Ns(r_nodes.size());

    GetGeometry().ShapeFunctionsValues(Ns, mHinges[i_hinge]);

    Vector r_coordinate(3);
    noalias(r_coordinate) = ZeroVector(3);
    for(std::size_t i = 0; i < r_nodes.size(); i++ )
        noalias(r_coordinate) += Ns[i]*r_nodes[i];

    return OversetCondition::PointType{r_coordinate[0], r_coordinate[1], r_coordinate[2]};
}

const Element & OversetCondition::rAdjacentElement() const
{ return * mpAdjacentElement; }

void OversetCondition::SetAdjacentElementAndSide( const Element * const p_adjacent_element, const std::size_t element_side )
{
    mpAdjacentElement = p_adjacent_element;
    mAdjacentElementSide = element_side;
}

const std::size_t OversetCondition::MeshBlockId() const
{
    return mpAdjacentElement->GetValue(BLOCK_ID);
}

}//namespace Kratos 