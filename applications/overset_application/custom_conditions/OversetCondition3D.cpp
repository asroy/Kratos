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
#include "custom_conditions/OversetCondition3D.h"

namespace Kratos
{
//default protected constructor
OversetCondition3D::OversetCondition3D(IndexType NewId)
    :   Condition(NewId),
        mIntegrationMethod{GetIntegrationMethod()},
        mAdjacentElementSide{0}
{}

OversetCondition3D::OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    :   Condition(NewId, pGeometry),
        mIntegrationMethod{GetIntegrationMethod()},
        mAdjacentElementSide{0}
{}

OversetCondition3D::OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :   Condition(NewId, pGeometry, pProperties),
        mIntegrationMethod{GetIntegrationMethod()},
        mAdjacentElementSide{0}
{}

OversetCondition3D::~OversetCondition3D()
{}

Condition::Pointer OversetCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new OversetCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void OversetCondition3D::GenerateHinges()
{
    const IntegrationPointsArrayType & r_integration_points = GetGeometry().IntegrationPoints(mIntegrationMethod);

    mHinge3Ds.clear();
    mHinge3Ds.reserve(r_integration_points.size());

    for( const IntegrationPointType & r_integration_point : r_integration_points )
        mHinge3Ds.push_back(Hinge3D{r_integration_point});
}

const Vector OversetCondition3D::HingeGlobalCoordinate(const std::size_t local_hinge_id ) const
{
    const GeometryType & r_nodes = GetGeometry();

    Vector Ns(r_nodes.size());

    GetGeometry().ShapeFunctionsValues(Ns, mHinge3Ds[local_hinge_id]);

    Vector r_coordinate(3);
    noalias(r_coordinate) = ZeroVector(3);
    for(std::size_t i = 0; i < r_nodes.size(); i++ )
        noalias(r_coordinate) += Ns[i]*r_nodes[i];

    return r_coordinate;
}

const std::vector<Hinge3D> & OversetCondition3D::Hinge3Ds() const
{ return mHinge3Ds; }

const Element::Pointer OversetCondition3D::pAdjacentElement() const
{ return mpAdjacentElement; }

void OversetCondition3D::SetAdjacentElementAndSide( const Element::WeakPointer & rp_adjacent_element, const std::size_t element_side )
{
    mpAdjacentElement = rp_adjacent_element.lock();
    mAdjacentElementSide = element_side;
}

const std::size_t OversetCondition3D::MeshBlockId() const
{
    return mpAdjacentElement->GetValue(BLOCK_ID);
}

Hinge3D & OversetCondition3D::rHinge3D(const std::size_t i)
{ return mHinge3Ds[i]; }

}//namespace Kratos 