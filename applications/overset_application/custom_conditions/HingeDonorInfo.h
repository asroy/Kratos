#if !defined(KRATOS_OVERSET_HINGEDONORINFO_H_INCLUDED )
#define  KRATOS_OVERSET_HINGEDONORINFO_H_INCLUDED

#include "includes/model_part.h"

namespace Kratos
{
namespace OversetAssembly
{

struct HingeDonorInfo
{
    std::size_t mModelPartId;
    ModelPart::ElementType::IndexType mElementId;
    std::vector<ModelPart::NodeType::IndexType> mNodesId;//mNodesId is for debugging
    Point<3> mBarycentricCoordinate;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
