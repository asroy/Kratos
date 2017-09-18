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
#include "includes/condition.h"
#include "includes/variables.h"

namespace Kratos
{

template<typename THinge>
class OversetCondition : public Condition
{
public:
    using Self = OversetCondition<THinge>;
    
    /// Counted pointer of OversetCondition
    KRATOS_CLASS_POINTER_DEFINITION(Self);

public:
    /// Default constructor
    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        :   Condition(NewId, pGeometry)
    {}

    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        :   Condition(NewId, pGeometry, pProperties)
    {}

    ~OversetCondition() override
    {}

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override
    {
        return Condition::Pointer new Self(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    void GenerateHinges();
    {
        Geometry::IntegrationPointsArrayType integration_points = pGetGeomrtry()->IntegrationPoints();

        mHinges.clear();
        mHinges.reserve(integration_points.size());

        for( const Geometry::IntegrationPointType & : r_integration_point )
            mHinges.push_back(THinge{r_integration_point});
    }

    const Element::Pointer pAdjacentElement()
    { return mpAdjacentElement; }

//member
private:
    std::vector<THinge> mHinges;
    Element::Pointer mpAdjacentElement;

friend class Serializer;

  	// A private default constructor necessary for serialization  
    OversetCondition()
        :   Condition()
    {}

};

} //namespace kratos 
#endif