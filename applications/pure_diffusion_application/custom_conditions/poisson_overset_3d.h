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

namespace Kratos
{
class PoissonOverset3D : public Condition
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Counted pointer of PointSource3D
    KRATOS_CLASS_POINTER_DEFINITION(PoissonOverset3D);
    
    /// Default constructor. 
    PoissonOverset3D(IndexType NewId, GeometryType::Pointer pGeometry);

    PoissonOverset3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    virtual ~PoissonOverset3D();

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);

protected:
 
 
private:
	// A private default constructor necessary for serialization  
    PoissonOverset3D() 
    :   Condition()
    {}

friend class Serializer;

};

} //namespace kratos 
#endif
