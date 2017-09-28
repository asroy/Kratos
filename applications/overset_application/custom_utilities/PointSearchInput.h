#if !defined(KRATOS_OVERSET_POINTSEARCHINPUT_H_INCLUDED )
#define  KRATOS_OVERSET_POINTSEARCHINPUT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct PointSearchInput
{
public:
    double mCoordinate[3];
    std::size_t mMeshBlockId;

private:
    void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mCoordinate[0]);
        r_serializer.Save(mCoordinate[1]);
        r_serializer.Save(mCoordinate[2]);
        r_serializer.Save(mMeshBlockId);
    }

    void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mCoordinate[0]);
        r_serializer.Load(mCoordinate[1]);
        r_serializer.Load(mCoordinate[2]);
        r_serializer.Load(mMeshBlockId);
    }

    void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{PointSearchInput: ";
        r_printer.Print(mCoordinate[0]);
        r_printer.Print(mCoordinate[1]);
        r_printer.Print(mCoordinate[2]);
        r_printer.Print(mMeshBlockId);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif