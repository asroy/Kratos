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
#include "includes/element.h"

#include "Hinge3D.h"

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
    OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

    OversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~OversetCondition3D() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    void GenerateHinges();

    const std::vector<Hinge3D> & Hinge3Ds() const;

    const std::vector<Vector> HingesGlobalCoordinate() const;

    const Element::WeakPointer pAdjacentElement() const;

    void SetAdjacentElementAndSide( const Element::WeakPointer & rp_adjacent_element, const std::size_t element_side );

protected:
    //default constructor necessary for serialization  
    OversetCondition3D(IndexType NewId = 0);

private:
    IntegrationMethod mIntegrationMethod;
    std::vector<Hinge3D> mHinge3Ds;
    Element::WeakPointer mpAdjacentElement;
    std::size_t mAdjacentElementSide;

friend class Serializer;
};

}//namespace Kratos 
#endif