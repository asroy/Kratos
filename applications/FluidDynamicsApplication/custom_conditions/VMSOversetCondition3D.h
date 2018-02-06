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

#if !defined(KRATOS_PoissonOverset3D_CONDITION_H_INCLUDED )
#define  KRATOS_PoissonOverset3D_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "overset_application/custom_conditions/OversetCondition.h"


namespace Kratos
{
class VMSOversetCondition3D : public OversetCondition
{
public:
    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(VMSOversetCondition3D);

    VMSOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

    VMSOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~VMSOversetCondition3D() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType & rLeftHandSideMatrix, VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo) override;
    
    void CalculateRightHandSide(VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo) override;

    void LocalEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo) override;

    void DonorEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;
 
protected:
	//default constructor necessary for serialization  
    VMSOversetCondition3D() 
        :   OversetCondition()
    {}

private:
    
friend class Serializer;

};

}//namespace Kratos 
#endif
