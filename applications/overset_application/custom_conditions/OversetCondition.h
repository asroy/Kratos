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

#if !defined(KRATOS_OVERSET_CONDITION_H_INCLUDED )
#define  KRATOS_OVERSET_CONDITION_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "includes/element.h"

#include "HingeDonorInfo.h"
#include "HingeDonorData.h"

namespace Kratos
{

class OversetCondition : public Condition
{
public:
    using PointType = Point;
    using IntegrationPointType = GeometryData::IntegrationPointType;
    using IntegrationPointsArrayType = GeometryData::IntegrationPointsArrayType;

    /// Counted pointer of OversetCondition
    KRATOS_CLASS_POINTER_DEFINITION(OversetCondition);

public:
    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    OversetCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~OversetCondition() override;

    Condition::Pointer Create(IndexType NewId, NodesArrayType const & ThisNodes, PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo) override;

    virtual void LocalEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo);

    virtual void DonorEquationIdVector(EquationIdVectorType & rResult, ProcessInfo & rCurrentProcessInfo);

    void GenerateHinges();

    const std::size_t NumberOfHinges() const;

    const IntegrationPointType & rHinge(const std::size_t i_hinge) const;

    OversetAssembly::HingeDonorInfo & rHingeDonorInfo(const std::size_t i_hinge);
    
    OversetAssembly::HingeDonorData & rHingeDonorData(const std::size_t i_hinge);

    const PointType HingeGlobalCoordinate(const std::size_t i_hinge) const;

    const Element & rAdjacentElement() const;

    void SetAdjacentElementAndSide( const Element * const p_adjacent_element, const std::size_t element_side );

    const std::size_t MeshBlockId() const;

    Vector ConditionJacobianToOutwardNormalVector(const Matrix & r_jacobian) const;//normal vector pointing outward
    
protected:
    //default constructor necessary for serialization  
    OversetCondition(IndexType NewId = 0);

private:
    const Element * mpAdjacentElement;
    std::size_t mAdjacentElementSide;

    std::vector<IntegrationPointType> mHinges;
    std::vector<OversetAssembly::HingeDonorInfo> mHingesDonorInfo;
    std::vector<OversetAssembly::HingeDonorData> mHingesDonorData;

friend class Serializer;
};

}//namespace Kratos 
#endif
