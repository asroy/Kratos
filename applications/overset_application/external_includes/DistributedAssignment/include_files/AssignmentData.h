#pragma once

namespace DistributedAssignment
{

template<typename TContractorKeyType,
         typename TAssignmentKeyType,
         typename TDataType>
class AssignmentData
{
public:
    AssignmentData()
    :   mAssignorKey(),
        mAssigneeKey(),
        mAssignmentKey(),
        mData()
    {}

    ~AssignmentData()
    {}

    TContractorKeyType GetAssignorKey() const
    { return mAssignorKey; }

    TContractorKeyType GetAssigneeKey() const
    { return mAssigneeKey; }

    TAssignmentKeyType GetAssignmentKey() const
    { return mAssignmentKey; }

    TDataType GetData() const
    { return mData; }

private:
    AssignmentData
    (   const TContractorKeyType & r_assignor_key,
        const TContractorKeyType & r_assignee_key,
        const TAssignmentKeyType & r_assignment_key,
        const TDataType & r_data )
    :   mAssignorKey(r_assignor_key),
        mAssigneeKey(r_assignee_key),
        mAssignmentKey(r_assignment_key),
        mData(r_data)
    {}

    void Save( DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mAssignorKey);
        r_serializer.Save(mAssigneeKey);
        r_serializer.Save(mAssignmentKey);
        r_serializer.Save(mData);
    }

    void Load( DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mAssignorKey);
        r_serializer.Load(mAssigneeKey);
        r_serializer.Load(mAssignmentKey);
        r_serializer.Load(mData);
    }

    void Profile( DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{AssignmentData: ";
        r_printer.Print(mAssignorKey);
        r_printer.Print(mAssigneeKey);
        r_printer.Print(mAssignmentKey);
        r_printer.Print(mData);
        std::cout << "},";
    }

    TContractorKeyType mAssignorKey;
    TContractorKeyType mAssigneeKey;
    TAssignmentKeyType mAssignmentKey;
    TDataType mData;


    template<   typename TAssignorType,
                typename TAssigneeType,
                typename TInputType,
                typename TOutputType,
                typename TCommunicatorType,
                template <typename TDummyType> class TDistributedContractorKeyIssuerType,
                template <typename TDummyType> class TDistributedAssignmentKeyIssuerType >
    friend class DistributedAssignmentManager;

    friend class DataUtility::Serializer;
    friend class DataUtility::DataProfile;
    friend class DataUtility::DataPrinter;
};

}