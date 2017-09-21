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

#if !defined(KRATOS_PointSource_3D_CONDITION_H_INCLUDED )
#define  KRATOS_PointSource_3D_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{
class PointSourceCondition3D
: public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PointSourceCondition3D);

    /// Default constructor. 
    PointSourceCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

    PointSourceCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~PointSourceCondition3D();

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;

private:
    // A private default constructor necessary for serialization  
    PointSourceCondition3D() 
    :   Condition()
    {}

friend class Serializer;

};

}//namespace Kratos 
#endif
