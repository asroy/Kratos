#if !defined(KRATOS_OVERSET_INTERPOLATIN_INPUT_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATIN_INPUT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

struct InterpolationInput
{
public:
    ModelPart::ElementType::IndexType mElementId;
    std::vector<ModelPart::NodeType::IndexType> mNodesId;
    double mBarycentricCoordinate[3];

private:
    void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mElementId);
        r_serializer.Save(mNodesId);
        r_serializer.Save(mBarycentricCoordinate[0]);
        r_serializer.Save(mBarycentricCoordinate[1]);
        r_serializer.Save(mBarycentricCoordinate[2]);
    }

    void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mElementId);
        r_serializer.Load(mNodesId);
        r_serializer.Load(mBarycentricCoordinate[0]);
        r_serializer.Load(mBarycentricCoordinate[1]);
        r_serializer.Load(mBarycentricCoordinate[2]);
    }

    void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{InterpolationInput: ";
        r_printer.Print(mElementId);
        r_printer.Print(mNodesId);
        r_printer.Print(mBarycentricCoordinate[0]);
        r_printer.Print(mBarycentricCoordinate[1]);
        r_printer.Print(mBarycentricCoordinate[2]);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
