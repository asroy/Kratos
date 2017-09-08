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

  class Poisson3D
	  : public Element
   {
   public:
     
     /// Counted pointer of Poisson3D
     KRATOS_CLASS_POINTER_DEFINITION(Poisson3D);


    /// Default constructor.
     Poisson3D(IndexType NewId, GeometryType::Pointer pGeometry);
     Poisson3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ Poisson3D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);


      
   protected:
   
   private:
	friend class Serializer;

       Poisson3D() : Element()
       {
       }
       
       
   }; // Class Poisson3D
}  // namespace Kratos.

#endif // KRATOS_POISSON_3D_ELEM_H_INCLUDED  defined
