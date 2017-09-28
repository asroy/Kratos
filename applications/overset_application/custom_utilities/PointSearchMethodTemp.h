#if !defined(KRATOS_OVERSET_POINT_SEARCH_METHOD_H_INCLUDED )
#define  KRATOS_OVERSET_POINT_SEARCH_METHOD_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "custom_utilities/DistributedAssignment.h"
#include "custom_utilities/SteSearcher.h"
#include "custom_utilities/PointSearchInput.h"
#include "custom_utilities/PointSearchOutput.h"

namespace Kratos
{
namespace OversetAssembly
{

template<typename TPointSearcher>
class PointSearchMethodTemp
{
public:
    using OversetCommunicator = DistributedAssignment::Communication::MpiCommunicator;
    using Location = OversetCommunicator::Location;
    using KeyIssuer = DistributedAssignment::DistributedAssignment::DistributedKeyIssuer<Location>;
    using Key = KeyIssuer::Key;

    using DummyModelPartHolder = DistributedAssignment::DistributedAssignment::DummyContractor<Key>;
    using DummyModelPartHolderManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<DummyModelPartHolder,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;
    
    using PointSearcherManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<TPointSearcher,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using PointSearchAssignmentManager = DistributedAssignment::DistributedAssignment::DistributedAssignmentManager<DummyModelPartHolder,TPointSearcher,PointSearchInput,PointSearchOutput,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using PointSearchAssignmentInputData  = typename PointSearchAssignmentManager::template AssignmentDataType<PointSearchInput>;
    using PointSearchAssignmentOutputData = typename PointSearchAssignmentManager::template AssignmentDataType<PointSearchOutput>;

    using BlockIdSet = std::set<std::size_t> ;
    using PointSearcherKeySetMapByMeshBlockId = std::map<std::size_t,std::set<Key,Key::LessThanComparator>>;

public:
    PointSearchMethodTemp() = delete;

    PointSearchMethodTemp
    ( 
        OversetCommunicator & r_communicator, 
        DummyModelPartHolderManager & r_dummy_model_part_holder_manager, 
        const std::size_t model_part_id,
        const ModelPart & r_model_part
    )
        :   mrOversetCommunicator{r_communicator},
            mrDummyModelPartHolderManager{r_dummy_model_part_holder_manager},
            mModelPartId{model_part_id},
            mpSearcherManager{nullptr},
            mpSearchAssignmentManager{nullptr}
    {
        //searcher
        std::vector<TPointSearcher *> local_point_searchers_pointer;
        TPointSearcher::BuildPointSearchersFromModelPart(r_model_part, mModelPartId, local_point_searchers_pointer);

        mpSearcherManager = new PointSearcherManager{mrOversetCommunicator};
        mpSearcherManager->RegisterLocalContractors( local_point_searchers_pointer, "SteSearcher" );
        mpSearcherManager->GenerateGlobalContractorsRegistry();

        //search assignment manager
        mpSearchAssignmentManager = new PointSearchAssignmentManager{r_communicator, mrDummyModelPartHolderManager, * mpSearcherManager};

        //block_id
        struct BlockIdAndSearcherKey
        {
        public:
            std::size_t mBlockId;
            Key mKey;

        private:
            void Save( DistributedAssignment::DataUtility::Serializer & r_serializer ) const
            {
                r_serializer.Save(mBlockId);
                r_serializer.Save(mKey);
            }
        
            void Load( DistributedAssignment::DataUtility::Serializer & r_serializer )
            {
                r_serializer.Load(mBlockId);
                r_serializer.Load(mKey);
            }
        
            void Profile( DistributedAssignment::DataUtility::DataProfile & r_profile ) const
            {
                r_profile.SetIsTrivial(false);
            }
        
            void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
            {
                std::cout << "{BlockIdAndSearcherKey: ";
                r_printer.Print(mBlockId);
                r_printer.Print(mKey);
                std::cout << "},";
            }
        
            friend class DistributedAssignment::DataUtility::Serializer;
            friend class DistributedAssignment::DataUtility::DataProfile;
            friend class DistributedAssignment::DataUtility::DataPrinter;
        };

        using BlockIdAndSearcherKeyVector = std::vector<BlockIdAndSearcherKey>;
        using BlockIdAndSearcherKeyVectorMapByLocation = OversetCommunicator::MapByLocationType<BlockIdAndSearcherKeyVector>;

        BlockIdAndSearcherKeyVector local_data_vector;
        BlockIdAndSearcherKeyVectorMapByLocation global_data_vector_map;

        for( const typename PointSearcherManager::ContractorPointerPairByContractorKey & r_contractor_pointer_pair : mpSearcherManager->LocalContractorsPointer() )
        {
            const Key key = r_contractor_pointer_pair.first;
            const TPointSearcher * const p_local_searcher = r_contractor_pointer_pair.second;
            const std::size_t block_id = p_local_searcher->MeshBlockId();
            local_data_vector.push_back({block_id, key});
        }

        mrOversetCommunicator.AllGather( local_data_vector, global_data_vector_map, 0 );

        for( const typename BlockIdAndSearcherKeyVectorMapByLocation::value_type & r_data_vector_pair : global_data_vector_map )
        {
            const BlockIdAndSearcherKeyVector & r_data_vector = r_data_vector_pair.second;

            for( const BlockIdAndSearcherKey & r_data : r_data_vector )
            {
                const std::size_t block_id = r_data.mBlockId;
                const Key key = r_data.mKey;

                mGlobalSearcherBlocksId.insert(block_id);
                mGlobalSearchersKeyForBlock[block_id].insert(key);
            }
        }
    }

    virtual ~PointSearchMethodTemp()
    {
        //searchers and manager
        for( typename PointSearcherManager::ContractorPointerPairByContractorKey & r_searcher_pair : mpSearcherManager->LocalContractorsPointer() )
            delete r_searcher_pair.second;

        delete mpSearcherManager;

        //assignment mananger
        delete mpSearchAssignmentManager;
    }

    const BlockIdSet & GlobalSearcherBlocksId() const
    { return mGlobalSearcherBlocksId; }

    const PointSearcherKeySetMapByMeshBlockId & GlobalSearchersKeyForBlock() const
    { return mGlobalSearchersKeyForBlock; }

    void ClearAllSearches()
    {
        mpSearchAssignmentManager->ClearAllAssignments();
    }

    Key AddSearch( const Key & r_searcher_key, const PointSearchInput & r_coordinate )
    {
        const Key & r_dummy_model_part_holder_key = *(mrDummyModelPartHolderManager.LocalContractorsKey().begin());

        return mpSearchAssignmentManager->AddAssignment( r_dummy_model_part_holder_key, r_searcher_key, r_coordinate );
    }

    void ExecuteAllSearches()
    {
        mpSearchAssignmentManager->ExecuteAllDistributedAssignments();
    }

    void GetSearchResults( std::vector<PointSearchAssignmentOutputData> & r_donors_info )
    {
        mpSearchAssignmentManager->GetResultsAtAssignor( r_donors_info );
    }

private:
    OversetCommunicator & mrOversetCommunicator;
    DummyModelPartHolderManager & mrDummyModelPartHolderManager;
    const std::size_t mModelPartId;
    PointSearcherManager * mpSearcherManager;
    PointSearchAssignmentManager * mpSearchAssignmentManager;
    BlockIdSet mGlobalSearcherBlocksId;
    PointSearcherKeySetMapByMeshBlockId mGlobalSearchersKeyForBlock;
};

}//namespace OversetAssemlby
}//namespace Kratos
#endif