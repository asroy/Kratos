#if !defined(KRATOS_OVERSET_INTERPOLATIONOUTPUT_H_INCLUDED )
#define  KRATOS_OVERSET_INTERPOLATIONOUTPUT_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{


class InterpolationOutput
{
public:
    using VariableKeyType = typename VariableData::KeyType;

public://make them public for debugging
    std::size_t mNumberOfDonorNodes;

    std::vector<double> mNs;
    std::vector<std::vector<double>> mDNsDXs;

    std::map<VariableKeyType,std::vector<std::size_t>> mDoubleVariablesNodesEquationId;
    std::map<VariableKeyType,double> mDoubleVariables;
    std::map<VariableKeyType,std::vector<double>> mDoubleVariablesDXs;

    std::map<VariableKeyType,std::vector<std::size_t>> mArray1dComponentVariablesNodesEquationId;
    std::map<VariableKeyType,double> mArray1dComponentVariables;
    std::map<VariableKeyType,std::vector<double>> mArray1dComponentVariablesDXs;

    std::vector<double> mCoordinate; //for debugging

public:
    InterpolationOutput()
        :   mNumberOfDonorNodes{0}
    {}

    virtual ~InterpolationOutput()
    {}

    void Clear()
    {
        mNumberOfDonorNodes = 0;
    
        mNs.clear();
        mDNsDXs.clear();

        mDoubleVariablesNodesEquationId.clear();
        mDoubleVariables.clear();
        mDoubleVariablesDXs.clear();

        mArray1dComponentVariablesNodesEquationId.clear();
        mArray1dComponentVariables.clear();
        mArray1dComponentVariablesDXs.clear();

        mCoordinate.clear();
    }

private:
    virtual void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mNumberOfDonorNodes);
        r_serializer.Save(mNs);
        r_serializer.Save(mDNsDXs);
        r_serializer.Save(mDoubleVariablesNodesEquationId);
        r_serializer.Save(mDoubleVariables);
        r_serializer.Save(mDoubleVariablesDXs);
        r_serializer.Save(mArray1dComponentVariablesNodesEquationId);
        r_serializer.Save(mArray1dComponentVariables);
        r_serializer.Save(mArray1dComponentVariablesDXs);
        r_serializer.Save(mCoordinate);
    }

    virtual void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mNumberOfDonorNodes);
        r_serializer.Load(mNs);
        r_serializer.Load(mDNsDXs);
        r_serializer.Load(mDoubleVariablesNodesEquationId);
        r_serializer.Load(mDoubleVariables);
        r_serializer.Load(mDoubleVariablesDXs);
        r_serializer.Load(mArray1dComponentVariablesNodesEquationId);
        r_serializer.Load(mArray1dComponentVariables);
        r_serializer.Load(mArray1dComponentVariablesDXs);
        r_serializer.Load(mCoordinate);
    }

    virtual void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    virtual void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{InterpolationOutput: ";
        r_printer.Print(mNumberOfDonorNodes);
        r_printer.Print(mNs);
        r_printer.Print(mDNsDXs);
        r_printer.Print(mDoubleVariablesNodesEquationId);
        r_printer.Print(mDoubleVariables);
        r_printer.Print(mDoubleVariablesDXs);
        r_printer.Print(mArray1dComponentVariablesNodesEquationId);
        r_printer.Print(mArray1dComponentVariables);
        r_printer.Print(mArray1dComponentVariablesDXs);
        r_printer.Print(mCoordinate);
        std::cout << "},";
    }

    friend class DistributedAssignment::DataUtility::Serializer;
    friend class DistributedAssignment::DataUtility::DataProfile;
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
