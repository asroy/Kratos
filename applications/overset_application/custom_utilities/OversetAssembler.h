#if !defined(KRATOS_OVERSET_ASSEMBLY_H_INCLUDED )
#define  KRATOS_OVERSET_ASSEMBLY_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_conditions/OversetCondition3D.h"
#include "custom_utilities/PointSearchMethodTemp.h"
#include "custom_utilities/SteSearcher.h"


namespace Kratos
{
namespace OversetAssembly
{


class OversetAssembler
{

public:
    using OversetCommunicator = DistributedAssignment::Communication::MpiCommunicator;
    using PointSearchMethod = PointSearchMethodTemp<SteSearcher>;

    OversetAssembler() = delete;
    
    OversetAssembler(const ModelPart & r_model_part)
        :   mrModelPart{r_model_part},
            mOversetCommunicator(),
            mpPointSearchMethod{new PointSearchMethod{mOversetCommunicator,mrModelPart}},
            mOversetCondition3Ds()
    {
        GetOversetConditionsFromInputModelPart();
    }

    virtual ~OversetAssembler()
    {
        delete mpPointSearchMethod;
    }

    void GetOversetConditionsFromInputModelPart()
    {
        // find out default overset condition
        mOversetCondition3Ds.clear();

        //have to cast away const of ModelPart here, because Conditions() is non const;
        const ModelPart::ConditionsContainerType & r_conditions_pointer = const_cast<ModelPart &>(mrModelPart).Conditions();

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = r_conditions_pointer.ptr_begin(); it_p_condition != r_conditions_pointer.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());
            if( p_overset_condition )
                mOversetCondition3Ds.push_back(* it_p_condition);
        }

        std::cout<<__func__<<": size mOversetCondition3Ds: "<<mOversetCondition3Ds.size()<<std::endl;

        // Generate overset condition-to-element adjacency
        using NodeIdVector = std::vector<std::size_t>;

        struct LessThanComparator
        {
            bool operator() ( const NodeIdVector & r_a_vector, const NodeIdVector & r_b_vector ) const
            {
                if( r_a_vector.size() < r_b_vector.size() )
                    return true;
                else if( r_a_vector.size() > r_b_vector.size() )
                    return false;

                std::set<std::size_t> a_set;
                std::set<std::size_t> b_set;
    
                for( const std::size_t & a : r_a_vector )
                    a_set.insert(a);

                for( const std::size_t & b : r_b_vector )
                    b_set.insert(b);

                std::set<std::size_t>::const_iterator it_a = a_set.begin();
                std::set<std::size_t>::const_iterator it_b = b_set.begin();

                while( it_a != a_set.end() )
                {
                    std::size_t a_node_id = *it_a;
                    std::size_t b_node_id = *it_b;

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
            Element::WeakPointer mpElement;
            std::size_t mElementSide;
        };

        using ConditionToElement = std::map<NodeIdVector,ElementAndSide,LessThanComparator>;

        //loop over elements' faces,
        //find key={condition_nodes} in map,
        //if not exsit, add {key={condition_nodes}, VALUE=ElementAndSide{element *, side}} into the map,
        //if already exist, delete the exisitng entry from the map
        ConditionToElement face_to_element_map;

        const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(mrModelPart).Elements();
        
        std::size_t num_face = 0;

        for( ModelPart::ElementsContainerType::ptr_const_iterator it_p_element = r_elements_pointer.ptr_begin(); it_p_element != r_elements_pointer.ptr_end(); it_p_element = std::next(it_p_element) )
        {
            const Element::GeometryType::GeometriesArrayType faces = (* it_p_element)->GetGeometry().Faces();

            for( std::size_t i = 0; i < faces.size(); i++ )
            {
                const Element::GeometryType & r_nodes = faces[i];

                NodeIdVector nodes_id(r_nodes.size());

                for( std::size_t j = 0; j < r_nodes.size(); j++ )
                    nodes_id[j] = r_nodes[j].GetId();

                ConditionToElement::const_iterator it = face_to_element_map.find(nodes_id);

                if( it == face_to_element_map.end() )
                    face_to_element_map[nodes_id] = {(* it_p_element),i};
                else
                    face_to_element_map.erase(it);

                num_face++;

                std::cout<<__func__<<": nodes_id: "<<nodes_id<<std::endl;
            }
        }

        std::cout<<__func__<<": num_face: "<<num_face<<std::endl;
        std::cout<<__func__<<": size face_to_element_map: "<<face_to_element_map.size()<<std::endl;
        

        //loop over overset overset conditions
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            const Condition::GeometryType & r_nodes = (* it_p_condition)->GetGeometry();

            NodeIdVector nodes_id(r_nodes.size());

            for( std::size_t i = 0; i < r_nodes.size(); i++ )
                nodes_id[i] = r_nodes[i].GetId();

            ConditionToElement::const_iterator it = face_to_element_map.find(nodes_id);

            if( it == face_to_element_map.end() )
            {
                //throw error please
                std::cout<<__func__<<": wrong! not found"<<std::endl;
                exit(EXIT_FAILURE);
            }

            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            if ( ! p_overset_condition )
            {
                //throw error please
                std::cout<<__func__<<": wrong! not OversetCondition3D"<<std::endl;
                exit(EXIT_FAILURE);
            }

            p_overset_condition->SetAdjacentElementAndSide( it->second.mpElement, it->second.mElementSide );
        }
    }

    const ModelPart::ConditionsContainerType & OversetCondition3Ds() const
    { return mOversetCondition3Ds; }

    void GenerateHinges()
    {
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            if ( ! p_overset_condition )
            {
                //throw error please
                std::cout<<__func__<<"wrong!"<<std::endl;
                exit(EXIT_FAILURE);
            }

            p_overset_condition->GenerateHinges();
        }
        
    }

    // //overset connectivity
    // void GenerateHingeDonorRelation()
    // {
    //     for( auto & rp_overset_condition : mOversetCondition3Ds )
    //     {}
    // }

    // void GetHingesValues()
    // {}
    
private:
    const ModelPart & mrModelPart;
    OversetCommunicator mOversetCommunicator;
    PointSearchMethod * const mpPointSearchMethod;
    ModelPart::ConditionsContainerType mOversetCondition3Ds;
};

}//namespace OverserAssembly
}//namespace Kratos
#endif