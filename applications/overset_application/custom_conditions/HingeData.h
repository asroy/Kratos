#if !defined(KRATOS_OVERSET_HINGEDATA_H_INCLUDED )
#define  KRATOS_OVERSET_HINGEDATA_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct HingeData
{
public:
    double mCoordinate[3];

private:
    void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mCoordinate[0]);
        r_serializer.Save(mCoordinate[1]);
        r_serializer.Save(mCoordinate[2]);
    }

    void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mCoordinate[0]);
        r_serializer.Load(mCoordinate[1]);
        r_serializer.Load(mCoordinate[2]);
    }

    void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{HingeData: ";
        r_printer.Print(mCoordinate[0]);
        r_printer.Print(mCoordinate[1]);
        r_printer.Print(mCoordinate[2]);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
