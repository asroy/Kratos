
namespace OversetAssembly
{

template<typename TCommunicator,
         typename <typename TDummy0> TPointSearcherType,
         template <typename TDummy1> TDistributedKeyIssuerType>
class PointSearchMethod
{
private:
    using Location = typename TCommunicator::Location;
    using DistributedKeyIssuer = typename TDistributedKeyIssuerType<Location>;
    using Key = typename DistributedKeyIssuer::Key;
    using PointSearcher = TPointSearcherType<Key>;

    using PointSearcherManager = DistributedAssignment::DistributedContractorManager<PointSearcher, TCommunicator, TDistributedKeyIssuerType>;
    using PointSearchAssignmentManager = DistributedAssignment::DistributedAssignmentManager<>
    using PointSearcherPointer = TPointSearcherType *;
    using PointSearcherPointerVector = std::vector<TPointSearcherTypePointer>;
    using DummyAssignor = DistributedAssignment::DummyContractor<>

public:
    PointSearchMethod( const TCommunicator & r_communicator, const & r_model_part )
        :   mrOversetCommunicator{r_communicator},
            mpDummyAssignorManager{nullptr},
            mpPointSearcherMananger{nullptr},
            mpPointSearchAssignmentManager{nullptr}
    {
        


        PointSearcherPointerVector local_point_searchers_pointer = BuildPointSearchers(r_model_part);

        mpPointSearcherManager = new PointSearcherManager{mrOversetCommunicator};
        RegisterLocalContractors( local_point_searchers_pointer, "PointSearcher" );
        GenerateGlobalContractorsRegistry();

        mpPointSearchAssignmentManager = new 
    }

    ~PointSearchMethod()
    {
        //PointSearcher(s) are constructed by the constructor of PointSearchMethod,
        //so we need to destruct them here
        PointSearcherPointerVector local_point_searchers_pointer = LocalContractorsPointer();
        DestroyPointSearchers( local_point_searchers_pointer );
    }

    //point searchers
    PointSearcherPointerVector BuildPointSearchers() const;
    
    void DestroyPointSearchers( PointSearcherPointerVector );

private:
    TCommunicator & mrOversetCommunicator;
    PointSearcherManager * mpPointSearcherManager;
    PointSearcherAssignmentManager * mpPointSearchAssignmentManager;
};


//specification for SteSearcher 
template<typename TCommunicator,
         template <typename TDummy> TDistributedKeyIssuerType>
std::vector<SteSearcher *> PointSearchMethod<SteSearcher,TCommunicator,TDistributedKeyIssuerType>::BuildPointSearchers() const
{
    using SteSearcherPointer = SteSearcher *;

    std::vector<SteSearcherPointer> point_searchers_pointer;

    //generate cnn, crd
    int num_nodes = mrModelPart.GetMesh().NumOfNodes();
    int num_elements = mrModelPart.GetMesh().NumOfElements();

    double * p_crd = new double [3*num_nodes];
    int * p_cnn = new int [4*num_elements];
    std::array<int> node_local_to_equation_id(num_nodes);
    std::map<int,int> node_equation_to_local_id(num_nodes);

    {
        int i = 0;
        for( auto & r_node : mrModelPart.GetMesh().Nodes() )
        {
            p_crd[3*i]   = r_node.X();
            p_crd[3*i+1] = r_node.Y();
            p_crd[3*i+2] = r_node.Z();
            i++;

            int equation_id = r_node.GetEquationId();

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

    //
    SteSearcherPointer p_point_searcher = new SteSearcher( p_crd, p_ccn, num_node, num_element );
    point_searchers_pointer.push_back(p_point_searcher);

    return point_searchers_pointer;
}

template<typename TCommunicator,
         template <typename TDummy> TDistributedKeyIssuerType>
void PointSearchMethod<SteSearcher,TCommunicator,TDistributedKeyIssuerType>::DestroyPointSearchers(std::vector<SteSearcher *> ste_searchers_pointer)
{
    for( SteSearcher * p_ste_searcher : ste_searchers_pointer )
        delete p_ste_searcher;
}

}