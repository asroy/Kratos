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

#if !defined(KRATOS_POISSON_HEATFLUX_3D_CONDITION_H_INCLUDED )
#define  KRATOS_POISSON_HEATFLUX_3D_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
class PoissonHeatFluxCondition3D : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PoissonHeatFluxCondition3D);

    /// Default constructor. 
    PoissonHeatFluxCondition3D(IndexType NewId, GeometryType::Pointer pGeometry);

    PoissonHeatFluxCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~PoissonHeatFluxCondition3D();

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

	void CalculateLeftHandSide(MatrixType & rLeftHandSideMatrix, ProcessInfo & rCurrentProcessInfo) override;
    
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;

private:
    // A private default constructor necessary for serialization  
    PoissonHeatFluxCondition3D() 
    :   Condition()
    {}

friend class Serializer;

};

}//namespace Kratos 
#endif
