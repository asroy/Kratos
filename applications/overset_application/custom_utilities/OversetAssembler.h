
namespace OversetAssembly
{

template<typename TOversetCommunicator,
         template <typename TDummy0> TPointSearcherType,
         template <typename TDummy1> TDistributedKeyIssuerType>
class OversetAssembler
{
private:
    ModelPart & mrModelPart;
    TOversetCommunicator mOversetCommunicator;
    OversetConditionContainer mOversetConditions;
    PointSearchMethod * mpPointSearchMethod;

public:
    using PointSearchMethod = TPointSearchMethod<TOversetCommunicator,TPointSearcherType,TDistributedKeyIssuerType>;
    using OversetConditionContainer = PointerVectorSet<OversetCondition, IndexedObject> ;

    OversetAssembler() = delete;
    
    OversetAssembler(const ModelPart & r_model_part)
        :   mrModelPart{r_model_part},
            mOversetCommunicator(),
            mOversetConditions(),
            mpPointSearchMethod{nullptr}
    {
        GetOversetConditionsFromInputModelPart();
        
        mpPointSearchMethod = new PointSearchMethod{mOversetCommunicator,mrModelPart};
    }

    ~OversetAssembler()
    {
        delete mpPointSearchMethod;
    }

    void GetOversetConditionsFromInputModelPart()
    {
        // find out default overset condition
        mOversetConditions.clear();

        for( auto p_condition : mrModelPart.GetMesh().Conditions() )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> p_condition.get();
            if( p_overset_condition );
                mOversetConditions.push_back(p_condition);
        }

        // Generate overset condition-to-element adjacency
        using NodeIdVector = std::vector<int>;

        struct LessThanComparator
        {
            bool operator() ( const NodeIdVector & a_vector, const NodeIdVector & b_vector ) const
            {
                if( a_vector.size() < b_vector.size() )
                    return true;
                else if( a_vector.size() > b_vector.size() )
                    return false;

                std::set<int> a_set;
                std::set<int> b_set;
    
                for( auto & a : a_vector )
                    a_set.insert(a);

                for( auto & b : b_vector )
                    b_set.insert(b);

                auto it_a = a_set.begin();
                auto it_b = b_set.begin();

                while( it_a != a_set.end() )
                {
                    int a_node_id = *it_a;
                    int b_node_id = *it_b;

                    if( a_node_id < b_node_id )
                        return true;

                    it_a = std::next(it_a);
                    it_b = std::next(it_b);
                }

                return false;
            }
        };

        struct ElementAndSide
        {
            Element::Pointer mpElement;
            int mElementSide;
        };

        using ConditionToElement = std::map<NodeIdVector,ElementAndSide,LessThanComparator>;

        //loop over elements' faces,
        //find key={condition_nodes} in map,
        //if not exsit, add {key={condition_nodes}, VALUE=ElementAndSide{element *, side}} into the map,
        //if already exist, delete the exisitng entry from the map
        ConditionToElement condition_to_element_map;

        auto & r_elements_pointer = mrModelPart.GetMesh().Elements();

        for( auto & rp_element : r_elements_pointer )
        {
            conditions = rp_element->conditions();

            for( int i = 0; it < conditions.size(); i++ )
            {
                nodes = conditions[i];

                NodeIdVector nodes_id(nodes.size());

                for( int j = 0; j < nodes.size(); j++ )
                    nodes_id[j] = nodes[j].GetEquationId();

                auto it = condition_to_element_map.find(nodes_id);

                if( it == condition_to_element_map.end() )
                    condition_to_element_map[nodes_id] = {rp_element,i};
                else
                    condition_to_element_map.erase(it);
            }
        }


        //loop over overset conditions, search for condition_nodes
        for( auto & rp_overset_condition : mOversetCondition )
        {
            auto & nodes = rp_overset_condition->Nodes();

            NodeIdVector nodes_id(nodes.size());

            for( int i = 0; i < nodes.size(); i++ )
                nodes_id[i] = nodes[i].GetEquationId();

            auto it = condition_to_element_map.find(nodes_id);

            if( it == condition_to_element_map.end() )
            {
                //throw error please
            }

            rp_overset_condition->mpElement = it->second.mpElement;
            rp_overset_condition->mElementSide = it->second.mElementSide;
        }
    }

    const OversetConditionContainer & OversetConditions() const
    { return mOversetConditions; }

    //
    void GenerateHinges()
    {
        for( auto & rp_overset_condition : mOversetConditions )
            rp_overset_condition->GenerateHinges();
    }

    //overset connectivity
    void GenerateHingeDonorRelation()
    {
        for( auto & rp_overset_condition : mOversetConditions )
        {
            for()
        }
    }

    //
    void GetHingesValues()
    {}

};

}