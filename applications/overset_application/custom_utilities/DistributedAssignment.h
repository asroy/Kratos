#if !defined(KRATOS_OVERSET_DISTRIBUTED_ASSIGNMENT_H_INCLUDED )
#define  KRATOS_OVERSET_DISTRIBUTED_ASSIGNMENT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{
namespace DistributedAssignment
{

#include "DistributedAssignment/Serializer.h"
#include "DistributedAssignment/DataProfile.h"
#include "DistributedAssignment/DataPrinter.h"
#include "DistributedAssignment/MpiCommunicator.h"
#include "DistributedAssignment/DistributedKeyIssuer.h"
#include "DistributedAssignment/DistributedContractorManager.h"
#include "DistributedAssignment/DistributedAssignmentManager.h"
#include "DistributedAssignment/AssignmentData.h"
#include "DistributedAssignment/DummyContractor.h"

}//namespace DistributedAssignment
}//namespace OversetAssembly
}//namespace Kratos
#endif
