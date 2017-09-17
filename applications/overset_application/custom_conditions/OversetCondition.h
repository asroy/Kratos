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
#include "Hinge.h"

namespace Kratos
{

template<std::szie_t TDimension,
         typename THingeData
         typename TData = double,
         typename TWeight = double>
class OversetCondition : public Condition
{
private:
    //type
    using Hinge = Hinge<TDimension, THingeData, TData, TWeight>;
    
    //member
    std::vector<Hinge> mHinges;
    Element::Pointer mpAdjacentElement; //smart pointer????

public:
    /// Counted pointer of OversetCondition
    KRATOS_CLASS_POINTER_DEFINITION(OversetCondition);

    /// Default constructor
    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~OversetCondition() override
    {}

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    const Element::Pointer pAdjacentElement()
    { return mpAdjacentElement; }

    void GenerateHinges();
    {
        mHinges.clear();

        int num_integration_point = //????

        mHinges.resize(num_of_integration_point);

        for( int i = 0; i < num_integration_point; i++ )
        {
            //assign data please
        }
    }

friend class Serializer;

  	// A private default constructor necessary for serialization  
    OversetCondition()
        :   Condition()
    {}

};

} //namespace kratos 
#endif