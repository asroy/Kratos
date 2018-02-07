#if !defined(KRATOS_OVERSET_HINGEDONORDATA_H_INCLUDED )
#define  KRATOS_OVERSET_HINGEDONORDATA_H_INCLUDED

namespace Kratos
{
namespace OversetAssembly
{

class HingeDonorData
{
public:
    using DoubleVariable = Variable<double>;
    using Array1dComponentVariable = VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>;
    using VariableKeyType = typename VariableData::KeyType;

public:
    bool mInitialized;

    std::size_t mNumberOfDonorNodes;

    std::vector<double> mNs;
    std::vector<std::vector<double>> mDNsDXs;

    std::map<VariableKeyType,std::vector<std::size_t>> mDoubleVariablesNodesEquationId;
    std::map<VariableKeyType,double> mDoubleVariables;
    std::map<VariableKeyType,std::vector<double>> mDoubleVariablesDXs;

    std::map<VariableKeyType,std::vector<std::size_t>> mArray1dComponentVariablesNodesEquationId;
    std::map<VariableKeyType,double> mArray1dComponentVariables;
    std::map<VariableKeyType,std::vector<double>> mArray1dComponentVariablesDXs;

    std::vector<double> mCoordinate;//for debugging

public:
    HingeDonorData()
        :   mInitialized{false},
            mNumberOfDonorNodes{0}
    {}

    ~HingeDonorData()
    {}

    void Clear()
    {
        mInitialized = false;

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

    bool IsInitialized() const
    { return mInitialized; }

    std::size_t NumberOfDonorNodes() const
    { return mNumberOfDonorNodes; }

    double GetN( const std::size_t i ) const
    { return mNs[i]; }

    std::vector<double> GetDNDXs( const std::size_t i ) const
    { return mDNsDXs[i]; }

    std::size_t GetDonorNodeEquationId( DoubleVariable & r_variable, const std::size_t i_node ) const
    { 
        auto it = mDoubleVariablesNodesEquationId.find( r_variable.Key() );

        if( it == mDoubleVariablesNodesEquationId.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mDoubleVariablesNodesEquationId"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return ( it->second )[i_node];
    }

    std::size_t GetDonorNodeEquationId( Array1dComponentVariable & r_variable, const std::size_t i_node ) const
    { 
        auto it = mArray1dComponentVariablesNodesEquationId.find( r_variable.Key() );

        if( it == mArray1dComponentVariablesNodesEquationId.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mArray1dComponentVariablesNodesEquationId"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return ( it->second )[i_node];
    }

    double GetValue( DoubleVariable & r_variable ) const
    { 
        auto it = mDoubleVariables.find( r_variable.Key() );

        if( it == mDoubleVariables.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mDoubleVariables"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }

    double GetValue( Array1dComponentVariable & r_variable ) const
    { 
        auto it = mArray1dComponentVariables.find( r_variable.Key() );

        if( it == mArray1dComponentVariables.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mArray1dComponentVariables"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }

    std::vector<double> GetDValueDXs( DoubleVariable & r_variable ) const
    { 
        auto it = mDoubleVariablesDXs.find( r_variable.Key() );

        if( it == mDoubleVariablesDXs.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mDoubleVariablesDXs"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }

    std::vector<double> GetDValueDXs( Array1dComponentVariable & r_variable ) const
    { 
        auto it = mArray1dComponentVariablesDXs.find( r_variable.Key() );

        if( it == mArray1dComponentVariablesDXs.end() )
        {
            std::cout<<__func__<<"wrong! key not find in mArray1dComponentVariablesDXs"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
