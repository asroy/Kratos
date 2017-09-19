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

#if !defined(KRATOS_OVERSET_CONDITION_H_INCLUDED )
#define  KRATOS_OVERSET_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "custom_conditions/Hinge3D.h"

namespace Kratos
{

class OversetCondition3D : public Condition
{
public:
    using IntegrationPointType = GeometryData::IntegrationPointType;
    using IntegrationPointsArrayType = GeometryData::IntegrationPointsArrayType;

    /// Counted pointer of OversetCondition3D
    KRATOS_CLASS_POINTER_DEFINITION(OversetCondition3D);

public:
    /// Default constructor
    OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
        :   Condition(NewId, pGeometry)
    {}

    OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        :   Condition(NewId, pGeometry, pProperties)
    {}

    ~OversetCondition3D() override
    {}

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer(new OversetCondition3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    void GenerateHinges()
    {
        IntegrationPointsArrayType integration_points = pGetGeometry()->IntegrationPoints();

        mHinges.clear();
        mHinges.reserve(integration_points.size());

        for( const IntegrationPointType & r_integration_point : integration_points )
            mHinges.push_back(Hinge3D{r_integration_point});
    }

    const Element::Pointer pAdjacentElement()
    { return mpAdjacentElement; }

//member
private:
    std::vector<Hinge3D> mHinges;
    Element::Pointer mpAdjacentElement;

  	// A private default constructor necessary for serialization  
    OversetCondition3D()
        :   Condition()
    {}

friend class Serializer;
};

}//namespace Kratos 
#endif