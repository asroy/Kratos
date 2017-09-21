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
#if !defined(KRATOS_POISSON_3D_ELEM_H_INCLUDED)
#define  KRATOS_POISSON_3D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

class PoissonElement3D
: public Element
{
public:
    /// Counted pointer of PoissonElement3D
    KRATOS_CLASS_POINTER_DEFINITION(PoissonElement3D);

    /// Default constructor.
    PoissonElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
    PoissonElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    ~PoissonElement3D() override;

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

protected:

private:
  	// A private default constructor necessary for serialization  
    PoissonElement3D() : Element()
    {}

friend class Serializer;

}; // Class PoissonElement3D
}  // namespace Kratos.

#endif // KRATOS_POISSON_3D_ELEM_H_INCLUDED  defined
