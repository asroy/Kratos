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

#if !defined(KRATOS_OVERSET_Hinge3D_H_INCLUDED )
#define  KRATOS_OVERSET_Hinge3D_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"

#include "integration/integration_point.h"

namespace Kratos
{

class Hinge3D : public IntegrationPoint<3>
{
private:
    using IntegrationPointType = IntegrationPoint<3>;

public:
    Hinge3D( const IntegrationPointType & r_integration_point )
        :   IntegrationPointType{r_integration_point}
    {}

    ~Hinge3D() override
    {}

private:
  	// A private default constructor necessary for serialization  
    Hinge3D()
        :   IntegrationPointType()
    {}

//member
private:
    std::vector<std::size_t> mDonorNodesEquationIds;

friend class Serializer;
};

}//namespace Kratos
#endif
