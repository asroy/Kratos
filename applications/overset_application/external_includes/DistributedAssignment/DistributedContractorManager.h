#pragma once
#include<iostream>
#include<vector>
#include<map>
#include<set>

namespace DistributedAssignment
{

template<   typename TContractorType,
            typename TCommunicatorType,
            template <typename TDummyType> class TDistributedContractorKeyIssuerType >
class DistributedContractorManager
{
public:
    using Location = typename TCommunicatorType::Location ;
    using ContractorKeyIssuer = TDistributedContractorKeyIssuerType<Location> ;
    using ContractorKey = typename ContractorKeyIssuer::Key ;
    using ContractorKeyVector = std::vector<ContractorKey> ;
    using ContractorKeySet = std::set<ContractorKey, typename ContractorKey::LessThanComparator> ;
    using ContractorPointer = TContractorType * ;
    using ContractorPointerMapByContractorKey = std::map<ContractorKey, ContractorPointer, typename ContractorKey::LessThanComparator> ;
    using LocationMapByContractorKey = std::map<ContractorKey, Location, typename ContractorKey::LessThanComparator> ;

    DistributedContractorManager() = delete;

    DistributedContractorManager( TCommunicatorType & communicator )
    :   mpCommunicator{& communicator},
        mContractorKeyIssuer(),
        mLocalContractorsKey(),
        mGlobalContractorsKey(),
        mLocalContractorsPointer(),
        mGlobalContractorsLocation()
    {}

    virtual ~DistributedContractorManager()
    {}

    void ClearContractorsRegistry()
    {
        mLocalContractorsKey.clear();
        mGlobalContractorsKey.clear();
        mLocalContractorsPointer.clear();
        mGlobalContractorsLocation.clear();
    }

    void RegisterLocalContractors( const std::vector<TContractorType *> & contractor_pointers, const std::string & r_name )
    {
        for( const TContractorType * const p_contractor : contractor_pointers )
            RegisterLocalContractor(*p_contractor, r_name);
    }

    void RegisterLocalContractor( TContractorType & r_contractor, const std::string & r_name )
    {
        const ContractorKey old_key = r_contractor.Key();

        if( old_key != ContractorKey::NoKey() )
        {
            DataUtility::DataPrinter printer;

            std::cout << __func__ << ": contractor ";

            printer.Print(old_key);

            std::cout <<"already has key! Do nothing" << std::endl;
        }

        const typename ContractorPointerMapByContractorKey::const_iterator it = mLocalContractorsPointer.find(old_key);

        if ( it == mLocalContractorsPointer.end() )
        {
            ContractorKey new_key = mContractorKeyIssuer.IssueNewKey();
            r_contractor.SetKey(new_key);
            r_contractor.SetName(r_name);
            mLocalContractorsPointer[new_key] = & r_contractor;
            mLocalContractorsKey.insert(new_key);
        }
        else
        {
            DataUtility::DataPrinter printer;

            std::cout << __func__ << ": contractor ";

            printer.Print(old_key);

            std::cout << "has already been registered locally! Do nothing" << std::endl;
        }
    }

    void GenerateGlobalContractorsRegistry( const int mpi_tag = 0 )
    {
        using ContractorKeyVectorMapByLocation = std::map<Location, ContractorKeyVector, typename Location::LessThanComparator> ;

        ContractorKeyVector local_contractor_key_vector;
        ContractorKeyVectorMapByLocation global_contractor_key_vector_map;

        //local contractors
        for( std::pair<const ContractorKey, ContractorPointer> & r_local_contractor_pointer_pair : mLocalContractorsPointer )
        {
            ContractorKey local_contractor_key = r_local_contractor_pointer_pair.first;
            local_contractor_key_vector.push_back(local_contractor_key);
        }

        //all gather
        mpCommunicator->AllGather( local_contractor_key_vector, global_contractor_key_vector_map, mpi_tag );

        //global contractor key to location map
        //global contractor key vector
        mGlobalContractorsLocation.clear();
        mGlobalContractorsKey.clear();

        for( std::pair<const Location, ContractorKeyVector> & r_global_contractor_key_vector_pair : global_contractor_key_vector_map )
        {
            Location location = r_global_contractor_key_vector_pair.first;
            ContractorKeyVector r_global_contractor_key_vector = r_global_contractor_key_vector_pair.second;

            for ( const ContractorKey & r_global_contractor_key : r_global_contractor_key_vector )
            {
                mGlobalContractorsLocation[r_global_contractor_key] = location;
                mGlobalContractorsKey.insert(r_global_contractor_key);
            }
        }
    }

    void PrintAllContractors() const
    {
        std::cout << __func__ << ": local contractors: " << std::endl;

        for( const std::pair<const ContractorKey, ContractorPointer> & r_contractor_pointer_pair : mLocalContractorsPointer )
        {
            const ContractorKey & r_contractor_key = r_contractor_pointer_pair.first;
            const TContractorType & r_contractor = *(r_contractor_pointer_pair.second);

            DataUtility::DataPrinter printer;
            printer.Print(r_contractor_key);
            printer.Print(r_contractor);

            std::cout<<std::endl;
        }

        std::cout << __func__ << ": global contractors: " << std::endl;

        for( const std::pair<const ContractorKey, Location> & r_contractor_location_pair : mGlobalContractorsLocation )
        {
            const ContractorKey & r_contractor_key = r_contractor_location_pair.first;
            const Location & r_location = r_contractor_location_pair.second;

            DataUtility::DataPrinter printer;
            printer.Print(r_contractor_key);
            printer.Print(r_location);

            std::cout<<std::endl;
        }
    }

    const ContractorKeySet & LocalContractorsKey() const
    { return mLocalContractorsKey; }

    const ContractorKeySet & GlobalContractorsKey() const
    { return mGlobalContractorsKey; }

    const ContractorPointerMapByContractorKey & LocalContractorsPointer() const
    { return mLocalContractorsPointer; }

    const LocationMapByContractorKey & GlobalContractorsLocation() const
    { return mGlobalContractorsLocation; }

    ContractorPointer const FindLocalContractorPointer( const ContractorKey key ) const
    {
        typename ContractorPointerMapByContractorKey::const_iterator it = mLocalContractorsPointer.find(key);
        if( it == mLocalContractorsPointer.end() )
            return nullptr;
        else
            return it->second;
    }

    Location FindGlobalContractorLocation( const ContractorKey key ) const
    {
        typename LocationMapByContractorKey::const_iterator it = mGlobalContractorsLocation.find(key);
        if( it == mGlobalContractorsLocation.end() )
            return Location::NoWhere();
        else
            return it->second;
    }

private:
    TCommunicatorType * const mpCommunicator;
    ContractorKeyIssuer mContractorKeyIssuer;

    ContractorKeySet mLocalContractorsKey;
    ContractorKeySet mGlobalContractorsKey;

    ContractorPointerMapByContractorKey mLocalContractorsPointer;
    LocationMapByContractorKey mGlobalContractorsLocation;
};

}
