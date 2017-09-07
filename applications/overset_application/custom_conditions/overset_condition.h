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

#if !defined(KRATOS_Overset_CONDITION_H_INCLUDED )
#define  KRATOS_Overset_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/variables.h"

namespace Kratos
{

template<std::size TDimension, typename TDataType = double, class TWeightType = double>
class Hinge : IntegrationPoint<TDimension, TDataType, TWeightType>
{
public:
    virtual ~Hinge();

protected:

private:
  	// A private default constructor necessary for serialization  
    Hinge();
  
friend class Serializer;
};


class OversetCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointSource2D
    KRATOS_CLASS_POINTER_DEFINITION(OversetCondition);

    /// Default constructor. 
    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    virtual ~OversetCondition();

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    Element & GetAdjacentElement();

protected:
    std::vector<Hinge> mHinges;
 
private:

friend class Serializer;

  	// A private default constructor necessary for serialization  
    OversetCondition()
        :   Condition()
    {}

};

} //namespace kratos 
#endif