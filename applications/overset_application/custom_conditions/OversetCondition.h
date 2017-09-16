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
#include "hinge.h"

namespace Kratos
{

template<std::szie_t TDimension,
         typename THingeDataType
         typename TDataType = double,
         typename TWeightType = double>
class OversetCondition : public Condition
{
public:
    using Hinge = Hinge<TDimension, THingeDataType, TDataType, TWeightType>;
    
    /// Counted pointer of OversetCondition
    KRATOS_CLASS_POINTER_DEFINITION(OversetCondition);

    /// Default constructor
    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~OversetCondition() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    const Element * AdjacentElementPointer()
    { return mpElement; }

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

protected:
    std::vector<Hinge> mHinges;
    Element::Pointer mpElement; //smart pointer????

private:

friend class Serializer;

  	// A private default constructor necessary for serialization  
    OversetCondition()
        :   Condition()
    {}

};

} //namespace kratos 
#endif