#if !defined(KRATOS_OVERSET_DONORINFO_H_INCLUDED )
#define  KRATOS_OVERSET_DONORINFO_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct DonorInfo
{
    std::vector<std::size_t> mDonorNodesEquationId;
    double mBarycentricCoordinate[3];
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
