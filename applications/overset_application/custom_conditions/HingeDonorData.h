#if !defined(KRATOS_OVERSET_HINGEDONORDATA_H_INCLUDED )
#define  KRATOS_OVERSET_HINGEDONORDATA_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct HingeDonorData
{
    std::vector<std::size_t> mEquationsId;
    std::vector<double> mNs;
    std::vector<std::vector<double>> mDNsDXs;
    double mTemperature;
    std::vector<double> mTempGradient;
    std::vector<double> mCoordinate;//for debugging
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
