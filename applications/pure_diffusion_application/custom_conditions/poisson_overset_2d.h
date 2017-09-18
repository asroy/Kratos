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

#if !defined(KRATOS_PoissonOverset2D_CONDITION_H_INCLUDED )
#define  KRATOS_PoissonOverset2D_CONDITION_H_INCLUDED

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
class PoissonOverset2D : public Condition
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Counted pointer of PointSource2D
    KRATOS_CLASS_POINTER_DEFINITION(PoissonOverset2D);
    
    /// Default constructor. 
    PoissonOverset2D(IndexType NewId, GeometryType::Pointer pGeometry);

    PoissonOverset2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~PoissonOverset2D() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;

protected:

 
private:
	// A private default constructor necessary for serialization  
    PoissonOverset2D() 
    :   Condition()
    {}

friend class Serializer;

};

} //namespace kratos 
#endif
