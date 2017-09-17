
namespace OversetAssembly
{

template<typename TOversetCommunicator,
         template <typename TDummy0> TPointSearcherType,
         template <typename TDummy1> TDistributedKeyIssuerType>
class PointSearchMethod
{
private:
    using Location = typename TOversetCommunicator::Location;
    using DistributedKeyIssuer = TDistributedKeyIssuerType<Location>;
    using Key = typename DistributedKeyIssuer::Key;

    using DummyContractor = DistributedAssignment::DummyContractor<Key>;
    using DummyContractorManager = DistributedAssignment::DistributedContractorManager<DummyContractor>

    using PointSearcher = TPointSearcherType<Key>;
    using PointSearcherManager = DistributedAssignment::DistributedContractorManager<PointSearcher,TOversetCommunicator,TDistributedKeyIssuerType>;

    using PointSearchAssignmentManager = DistributedAssignment::DistributedAssignmentManager<DummyContractor,PointSearcher,Coordinate,DonorInfo,TOversetCommunicator,TDistributedKeyIssuerType,TDistributedKeyIssuerType>

public:
    PointSearchMethod( const TOversetCommunicator & r_communicator, const & r_model_part )
        :   mrOverseTOversetCommunicator{r_communicator},
            mpDummyAssignorManager{nullptr},
            mpPointSearcherManager{nullptr},
            mpPointSearchAssignmentManager{nullptr}
    {
        //dummy assignor
        mpDummyAssignorManager = new DummyContractorManager{mrOverseTOversetCommunicator};
        mpDummyAssignorManager->RegisterLocalContractor( new DummyContractor() );
        mpDummyAssignorManager->GenerateGlobalContractorsRegistry();

        //searcher
        std::vector<TPointSearcherType *> local_point_searchers_pointer = BuildPointSearchers(r_model_part);
        mpPointSearcherManager = new PointSearcherManager{mrOverseTOversetCommunicator};
        mpPointSearcherManager->RegisterLocalContractors( local_point_searchers_pointer, "PointSearcher" );
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
        for( PointSearcherPointer p_searcher : mpPointSearcherManager->LocalContractorsPointer() )
            delete p_searcher;

        delete mpPointSearcherManager;

        //assignment mananger
        delete mpPointSearchAssignmentManager;
    }
    
private:
    TOversetCommunicator & mrOverseTOversetCommunicator;
    DummyContractorManager * mpDummyAssignorManager;
    PointSearcherManager * mpPointSearcherManager;
    PointSearcherAssignmentManager * mpPointSearchAssignmentManager;
};

//specialization for SteSearcher
template<typename TOversetCommunicator,
         template <typename TDummy1> TDistributedKeyIssuerType>
std::vector<SteSearcher *> PointSearchMethod<TOversetCommunicator,SteSearcher,TDistributedKeyIssuerType>::BuildPointSearchers() const
{
    std::vector<SteSearcher *> point_searchers_pointer;

    //generate cnn, crd
    int num_node = mrModelPart.GetMesh().NumOfNodes();
    int num_element = mrModelPart.GetMesh().NumOfElements();

    double * p_crd = new double [3*num_node];
    int * p_cnn = new int [4*num_element];
    std::vector<std::size_t> node_local_to_equation_id(num_node);
    std::map<std::size_t,int> node_equation_to_local_id;

    {
        int i = 0;
        for( auto & r_node : mrModelPart.GetMesh().Nodes() )
        {
            p_crd[3*i]   = r_node.X();
            p_crd[3*i+1] = r_node.Y();
            p_crd[3*i+2] = r_node.Z();
            i++;

            std::size_t equation_id = r_node.GetEquationId();

            node_local_to_equation_id[i] = equation_id;
            node_equation_to_local_id[equation_id] = i;
        }
    }

    {
        int i = 0;
        for( auto & r_element : mrModelPart.GetMesh().Elements() )
        {
            assert( r_element.Nodes().size() == 4 )

            for( auto r_node : r_element.Nodes() )
            {
                p_cnn[i] = node_equation_to_local_id[r_node.GetEquationId()];
                i++;
            }
        }
    }
    
    //create SteSearcher
    SteSearcher * p_point_searcher = new SteSearcher( p_crd, p_ccn, node_local_to_equation_id, num_node, num_element );
    point_searchers_pointer.push_back(p_point_searcher);

    return point_searchers_pointer;
}


}