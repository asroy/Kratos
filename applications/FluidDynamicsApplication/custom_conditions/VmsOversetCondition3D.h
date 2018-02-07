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

#if !defined(KRATOS_VMS_OVERSET_CONDITION_3D_H_INCLUDED )
#define  KRATOS_VMS_OVERSET_CONDITION_3D_H_INCLUDED

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
class VmsOversetCondition3D : public OversetCondition
{
public:
    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(VmsOversetCondition3D);

    VmsOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

    VmsOversetCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~VmsOversetCondition3D() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType & rLeftHandSideMatrix, VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo) override;
    
    void CalculateRightHandSide(VectorType & rRightHandSideVector, ProcessInfo & rCurrentProcessInfo) override;

    void LocalEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo) override;

    void DonorEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;
 
protected:
	//default constructor necessary for serialization  
    VmsOversetCondition3D(IndexType NewId = 0);

private:
    
friend class Serializer;

};

}//namespace Kratos 
#endif
