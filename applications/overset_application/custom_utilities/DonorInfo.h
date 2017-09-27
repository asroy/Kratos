#if !defined(KRATOS_OVERSET_DONORINFO_H_INCLUDED )
#define  KRATOS_OVERSET_DONORINFO_H_INCLUDED

#include <custom_utilities/Coordinate.h>

namespace Kratos
{
namespace OversetAssembly
{

struct DonorInfo
{
public:
    bool mFound;
    std::size_t mDonorModelPartId;
    std::size_t mDonorElementId;
    std::vector<std::size_t> mDonorNodesId;
    double mDonorBarycentricCoordinate[3];

    double mDistance;
    std::size_t mDonorMeshBlockId;
    Coordinate mInterpolatedCoordinate;

private:
    void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mFound);
        r_serializer.Save(mDonorModelPartId);
        r_serializer.Save(mDonorElementId);
        r_serializer.Save(mDonorNodesId);
        r_serializer.Save(mDonorBarycentricCoordinate[0]);
        r_serializer.Save(mDonorBarycentricCoordinate[1]);
        r_serializer.Save(mDonorBarycentricCoordinate[2]);
        r_serializer.Save(mDistance);
        r_serializer.Save(mDonorMeshBlockId);
        r_serializer.Save(mInterpolatedCoordinate);
    }

    void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mFound);
        r_serializer.Load(mDonorModelPartId);
        r_serializer.Load(mDonorElementId);
        r_serializer.Load(mDonorNodesId);
        r_serializer.Load(mDonorBarycentricCoordinate[0]);
        r_serializer.Load(mDonorBarycentricCoordinate[1]);
        r_serializer.Load(mDonorBarycentricCoordinate[2]);
        r_serializer.Load(mDistance);
        r_serializer.Load(mDonorMeshBlockId);
        r_serializer.Load(mInterpolatedCoordinate);
    }

    void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{DonorInfo: ";
        r_printer.Print(mFound);
        r_printer.Print(mDonorModelPartId);
        r_printer.Print(mDonorElementId);
        r_printer.Print(mDonorNodesId);
        r_printer.Print(mDonorBarycentricCoordinate[0]);
        r_printer.Print(mDonorBarycentricCoordinate[1]);
        r_printer.Print(mDonorBarycentricCoordinate[2]);
        r_printer.Print(mDistance);
        r_printer.Print(mDonorMeshBlockId);
        r_printer.Print(mInterpolatedCoordinate);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
