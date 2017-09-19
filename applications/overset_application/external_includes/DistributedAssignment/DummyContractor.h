#pragma once
#include<string>
#include<set>

namespace DistributedAssignment
{

template<typename TContractorKeyType>
class DummyContractor
{
public:
    using ContractorKeySet = std::set<TContractorKeyType, typename TContractorKeyType::LessThanComparator>;

    DummyContractor()
    :   mName(),
        mKey{TContractorKeyType::NoKey()}
    {}

    ~DummyContractor()
    {}

    void SetName( const std::string & r_name )
    { mName = r_name; }

    std::string Name() const
    { return mName; }

    void SetKey(const TContractorKeyType contrator_key)
    { mKey = contrator_key; }

    TContractorKeyType Key() const
    { return mKey; }

    void Execute( const int & r_in, int & r_out )
    {
        r_out = r_in;
    }

private:
    void Print( const DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{DummyContractor: ";
        std::cout << "{Name: "<< mName <<"},",
        std::cout << "{ContractorKey: ";
        r_printer.Print(mKey);
        std::cout << "},";
        std::cout << "},";
    }

    std::string mName;
    TContractorKeyType mKey;

    friend class DataUtility::DataPrinter;
};

}