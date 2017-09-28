#if !defined(KRATOS_OVERSET_POINTSEARCHOUTPUT_H_INCLUDED )
#define  KRATOS_OVERSET_POINTSEARCHOUTPUT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct PointSearchOutput
{
public:
    bool mFound;
    std::size_t mModelPartId;
    std::size_t mElementId;
    std::vector<std::size_t> mNodesId;
    double mBarycentricCoordinate[3];

    double mDistance;
    std::size_t mMeshBlockId;
    double mInterpolatedCoordinate[3];

private:
    void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mFound);
        r_serializer.Save(mModelPartId);
        r_serializer.Save(mElementId);
        r_serializer.Save(mNodesId);
        r_serializer.Save(mBarycentricCoordinate[0]);
        r_serializer.Save(mBarycentricCoordinate[1]);
        r_serializer.Save(mBarycentricCoordinate[2]);
        r_serializer.Save(mDistance);
        r_serializer.Save(mMeshBlockId);
        r_serializer.Save(mInterpolatedCoordinate[0]);
        r_serializer.Save(mInterpolatedCoordinate[1]);
        r_serializer.Save(mInterpolatedCoordinate[2]);
    }

    void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mFound);
        r_serializer.Load(mModelPartId);
        r_serializer.Load(mElementId);
        r_serializer.Load(mNodesId);
        r_serializer.Load(mBarycentricCoordinate[0]);
        r_serializer.Load(mBarycentricCoordinate[1]);
        r_serializer.Load(mBarycentricCoordinate[2]);
        r_serializer.Load(mDistance);
        r_serializer.Load(mMeshBlockId);
        r_serializer.Load(mInterpolatedCoordinate[0]);
        r_serializer.Load(mInterpolatedCoordinate[1]);
        r_serializer.Load(mInterpolatedCoordinate[2]);
    }

    void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{PointSearchOutput: ";
        r_printer.Print(mFound);
        r_printer.Print(mModelPartId);
        r_printer.Print(mElementId);
        r_printer.Print(mNodesId);
        r_printer.Print(mBarycentricCoordinate[0]);
        r_printer.Print(mBarycentricCoordinate[1]);
        r_printer.Print(mBarycentricCoordinate[2]);
        r_printer.Print(mDistance);
        r_printer.Print(mMeshBlockId);
        r_printer.Print(mInterpolatedCoordinate[0]);
        r_printer.Print(mInterpolatedCoordinate[1]);
        r_printer.Print(mInterpolatedCoordinate[2]);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
