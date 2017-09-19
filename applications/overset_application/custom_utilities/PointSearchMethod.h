#if !defined(KRATOS_OVERSET_POINT_SEARCH_METHOD_H_INCLUDED )
#define  KRATOS_OVERSET_POINT_SEARCH_METHOD_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "DistributedAssignment.h"
#include "SteSearcher.h"
#include "Coordinate.h"
#include "DonorInfo.h"

namespace Kratos
{
namespace OversetAssembly
{

template<typename TPointSearcher>
class PointSearchMethod
{
private:
    using OversetCommunicator = DistributedAssignment::Communication::MpiCommunicator;
    using Location = typename OversetCommunicator::Location;
    using KeyIssuer = DistributedAssignment::DistributedAssignment::DistributedKeyIssuer<Location>;
    using Key = typename KeyIssuer::Key;

    using DummyContractor = DistributedAssignment::DistributedAssignment::DummyContractor<Key>;
    using DummyContractorManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<DummyContractor,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using PointSearcherManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<TPointSearcher,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using PointSearchAssignmentManager = DistributedAssignment::DistributedAssignment::DistributedAssignmentManager<DummyContractor,TPointSearcher,Coordinate,DonorInfo,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

public:
    PointSearchMethod() = delete;

    PointSearchMethod( const OversetCommunicator & r_communicator, const ModelPart & r_model_part )
        :   mrOverseTOversetCommunicator{r_communicator},
            mpDummyAssignorManager{nullptr},
            mpPointSearcherManager{nullptr},
            mpPointSearchAssignmentManager{nullptr}
    {
        //dummy assignor
        mpDummyAssignorManager = new DummyContractorManager{mrOverseTOversetCommunicator};
        mpDummyAssignorManager->RegisterLocalContractor( * (new DummyContractor()), "DummyContractor" );
        mpDummyAssignorManager->GenerateGlobalContractorsRegistry();

        //searcher
        std::vector<TPointSearcher *> local_point_searchers_pointer;
        TPointSearcher::BuildPointSearchersFromModelPart(mrModelPart, local_point_searchers_pointer);

        mpPointSearcherManager = new PointSearcherManager{mrOverseTOversetCommunicator};
        mpPointSearcherManager->RegisterLocalContractors( local_point_searchers_pointer, "SteSearcher" );
        mpPointSearcherManager->GenerateGlobalContractorsRegistry();

        //search assignment manager
        mpPointSearchAssignmentManager = new PointSearchAssignmentManager{r_communicator, *mpDummyAssignorManager, *mpPointSearcherManager};
    }

    ~PointSearchMethod()
    {
        //dummy assignor and manager
        for( DummyContractor * p_dummy_assignor : mpDummyAssignorManager->LocalContractorsPointer() )
            delete p_dummy_assignor;

        delete mpDummyAssignorManager;
        
        //searchers and manager
        for( TPointSearcher * p_searcher : mpPointSearcherManager->LocalContractorsPointer() )
            delete p_searcher;

        delete mpPointSearcherManager;

        //assignment mananger
        delete mpPointSearchAssignmentManager;
    }

    void ClearAllSearches()
    {
        mpPointSearchAssignmentManager->ClearAllAssignments();
    }

    void AddSearch( const Key & r_searcher_key, const Coordinate & r_coordinate )
    {
        mpPointSearchAssignmentManager->AddAssignment( Key::NoKey(), r_searcher_key, r_coordinate );
    }

    void ExecuteAllSearches()
    {
        mpPointSearchAssignmentManager->ExecuteAllDistributedAssignments();
    }

    void GetSearchResults( std::vector<DonorInfo> & r_donors_info )
    {
        mpPointSearchAssignmentManager->GetResultsAtAssignor( r_donors_info );
    }

private:
    OversetCommunicator & mrOverseTOversetCommunicator;
    DummyContractorManager * mpDummyAssignorManager;
    PointSearcherManager * mpPointSearcherManager;
    PointSearchAssignmentManager * mpPointSearchAssignmentManager;
};

}//namespace OversetAssemlby
}//namespace Kratos
#endif