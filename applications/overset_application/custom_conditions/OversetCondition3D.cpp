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
#include "custom_conditions/OversetCondition3D.h"

namespace Kratos
{
OversetCondition3D::OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    :   Condition(NewId, pGeometry)
{}

OversetCondition3D::OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    :   Condition(NewId, pGeometry, pProperties)
{}

OversetCondition3D::~OversetCondition3D()
{}

Condition::Pointer OversetCondition3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new OversetCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void OversetCondition3D::GenerateHinges()
{
    IntegrationPointsArrayType integration_points = pGetGeometry()->IntegrationPoints();

    mHinge3Ds.clear();
    mHinge3Ds.reserve(integration_points.size());

    for( const IntegrationPointType & r_integration_point : integration_points )
        mHinge3Ds.push_back(Hinge3D{r_integration_point});
}

const std::vector<Hinge3D> & OversetCondition3D::Hinge3Ds() const
{ return mHinge3Ds; }

const Element::WeakPointer OversetCondition3D::pAdjacentElement() const
{ return mpAdjacentElement; }

void OversetCondition3D::SetAdjacentElementAndSide( const Element::WeakPointer & rp_adjacent_element, const std::size_t element_side )
{
    mpAdjacentElement = rp_adjacent_element;
    mAdjacentElementSide = element_side;
}

}//namespace Kratos 