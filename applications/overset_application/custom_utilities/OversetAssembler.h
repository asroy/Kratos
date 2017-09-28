#if !defined(KRATOS_OVERSET_ASSEMBLY_H_INCLUDED )
#define  KRATOS_OVERSET_ASSEMBLY_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_conditions/OversetCondition.h"
#include "custom_utilities/DistributedAssignment.h"
#include "custom_utilities/PointSearchMethodTemp.h"
#include "custom_utilities/SteSearcher.h"
#include "custom_utilities/InterpolationMethod.h"
#include "custom_utilities/Interpolator.h"


namespace Kratos
{
namespace OversetAssembly
{


class OversetAssembler
{

private:
    using OversetCommunicator = DistributedAssignment::Communication::MpiCommunicator;
    using Location = OversetCommunicator::Location;
    using KeyIssuer = DistributedAssignment::DistributedAssignment::DistributedKeyIssuer<Location>;
    using Key = KeyIssuer::Key;
    using DummyModelPartHolder = DistributedAssignment::DistributedAssignment::DummyContractor<Key>;
    using DummyModelPartHolderManager = DistributedAssignment::DistributedAssignment::DistributedContractorManager<DummyModelPartHolder,OversetCommunicator,DistributedAssignment::DistributedAssignment::DistributedKeyIssuer>;

    using PointSearchMethod = PointSearchMethodTemp<SteSearcher>;
    using PointSearcherKey = PointSearchMethod::Key;
    using PointSearchAssignmentKey = PointSearchMethod::Key;
    
    using InterpolatorKey = InterpolationMethod::Key;
    using InterpolationAssignmentKey = InterpolationMethod::Key;

    using PointType = Point<3>;

public:
    OversetAssembler() = delete;

    OversetAssembler(const ModelPart & r_model_part)
        :   mrModelPart{r_model_part},
            mOversetCommunicator(),
            mDummyModelPartHolder(),
            mDummyModelPartHolderManager{mOversetCommunicator},
            mpPointSearchMethod{nullptr},
            mpInterpolationMethod{nullptr},
            mOversetConditions()
    {
        //dummy model part holder
        mDummyModelPartHolderManager.RegisterLocalContractor( mDummyModelPartHolder, "DummyModelPartHolder" );
        mDummyModelPartHolderManager.GenerateGlobalContractorsRegistry();

        {
            std::size_t model_part_id = 0;
            for ( const Key & r_model_part_holder_key : mDummyModelPartHolderManager.GlobalContractorsKey() )
            {
                mModelPartKeyToId[r_model_part_holder_key] = model_part_id;
                mModelPartIdToKey[model_part_id] = r_model_part_holder_key;
                model_part_id++;
            }
        }

        //point search method
        {
            if( mModelPartKeyToId.find(mDummyModelPartHolder.Key()) == mModelPartKeyToId.end() )
            {
                std::cout<<__func__<<"wrong! model part holder not found"<<std::endl;
                exit(EXIT_FAILURE);
            }

            const std::size_t model_part_id = mModelPartKeyToId[mDummyModelPartHolder.Key()];
            mpPointSearchMethod = new PointSearchMethod{mOversetCommunicator,mDummyModelPartHolderManager, model_part_id, mrModelPart};
        }

        //interpolation method
        mpInterpolationMethod = new InterpolationMethod( mOversetCommunicator, mDummyModelPartHolderManager, mrModelPart );
        
        //overset conditions
        GenerateOversetConditionsFromInputModelPart();
    }

    virtual ~OversetAssembler()
    {
        delete mpPointSearchMethod;
        delete mpInterpolationMethod;
    }

    void GenerateOversetConditionsFromInputModelPart()
    {
        // find out default overset condition
        mOversetConditions.clear();

        //have to cast away const of ModelPart here, because Conditions() is non const;
        const ModelPart::ConditionsContainerType & r_conditions_pointer = const_cast<ModelPart &>(mrModelPart).Conditions();

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = r_conditions_pointer.ptr_begin(); it_p_condition != r_conditions_pointer.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());
            if( p_overset_condition )
                mOversetConditions.push_back(* it_p_condition);
        }

        std::cout<<__func__<<": size mOversetConditions: "<<mOversetConditions.size()<<std::endl;

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
                    else if( a_node_id > b_node_id )
                        return false;

                    it_a = std::next(it_a);
                    it_b = std::next(it_b);
                }

                return false;
            }
        };

        struct ElementAndSide
        {
            const Element * mpElement;
            std::size_t mElementSide;
        };

        using ConditionToElement = std::map<NodeIdVector,ElementAndSide,LessThanComparator>;

        //loop over elements' faces,
        //find key={condition_nodes} in map,
        //if not exsit, add {key={condition_nodes}, VALUE=ElementAndSide{element *, side}} into the map,
        //if already exist, delete the exisitng entry from the map
        ConditionToElement face_to_element_map;

        const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(mrModelPart).Elements();
        
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
                    face_to_element_map[nodes_id] = {(* it_p_element).get(),i};
                else
                    face_to_element_map.erase(it);
            }
        }

        std::cout<<__func__<<": size face_to_element_map: "<<face_to_element_map.size()<<std::endl;


        //loop over overset conditions
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
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

            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            if ( ! p_overset_condition )
            {
                //throw error please
                std::cout<<__func__<<": wrong! not OversetCondition"<<std::endl;
                exit(EXIT_FAILURE);
            }

            p_overset_condition->SetAdjacentElementAndSide( it->second.mpElement, it->second.mElementSide );
        }
    }

    const ModelPart::ConditionsContainerType & rOversetConditions() const
    { return mOversetConditions; }

    void GenerateHinges()
    {
        std::size_t num_overset_condition = 0;
        std::size_t num_hinge = 0;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            if ( ! p_overset_condition )
            {
                //throw error please
                std::cout<<__func__<<"wrong! not OversetCondition *"<<std::endl;
                exit(EXIT_FAILURE);
            }

            p_overset_condition->GenerateHinges();

            num_overset_condition++;
            num_hinge += p_overset_condition->NumberOfHinges();
        }

        std::cout<<__func__<<": num_overset_condition "<<num_overset_condition<<", num_hinge "<<num_hinge<<std::endl;
    }

    void SearchHingesDonor()
    {
        //
        struct HingeKey
        {
            OversetCondition::IndexType mConditionId;
            std::size_t mHingLocalId;

            bool operator< (const HingeKey & r_other ) const
            {
                if( mConditionId < r_other.mConditionId )
                    return true;
                else if ( mConditionId > r_other.mConditionId )
                    return false;
                else if ( mHingLocalId < r_other.mHingLocalId )
                    return true;

                return false;
            }
        };

        using HingeToAssignmentVectorMap = std::map<HingeKey,std::vector<PointSearchAssignmentKey>>;

        HingeToAssignmentVectorMap hinge_to_assignments_map;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            OversetCondition::IndexType condition_id = p_overset_condition->GetId();

            //add search assignment
            std::size_t condition_block_id = p_overset_condition->MeshBlockId();

            for ( const PointSearchMethod::PointSearcherKeySetMapByMeshBlockId::value_type & r_pair : mpPointSearchMethod->GlobalSearchersKeyForBlock() )
            {
                std::size_t searcher_block_id = r_pair.first;

                if( condition_block_id != searcher_block_id )
                {
                    for( const PointSearcherKey & r_searcher_key : r_pair.second )
                    {
                        for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->NumberOfHinges(); i_hinge++ )
                        {  
                            const Vector hinge_coordinate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                            //add search assignment
                            PointSearchAssignmentKey search_assignment_key = mpPointSearchMethod->AddSearch( r_searcher_key, {hinge_coordinate[0], hinge_coordinate[1], hinge_coordinate[2]} );

                            //associate search assigment with hinge
                            HingeKey hinge_key{condition_id,i_hinge};
                            hinge_to_assignments_map[hinge_key].push_back(search_assignment_key);
                        }
                    }
                }
            }
        }

        //execute search
        mpPointSearchMethod->ExecuteAllSearches();

        //get search result
        std::vector<PointSearchMethod::PointSearchAssignmentOutputData> donor_info_data_vector;
        mpPointSearchMethod->GetSearchResults( donor_info_data_vector );

        {
            DistributedAssignment::DataUtility::DataPrinter printer;
            printer.Print(donor_info_data_vector);
        }

        //get search result mapped by assignment key
        std::map<PointSearchAssignmentKey,PointSearchOutput,PointSearchAssignmentKey::LessThanComparator> donor_info_map;

        for( const auto & r_donor_info_data : donor_info_data_vector )
        {
            donor_info_map[r_donor_info_data.GetAssignmentKey()] = r_donor_info_data.GetData();
        }

        // get donor result for hinges
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            OversetCondition::IndexType condition_id = p_overset_condition->GetId();

            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->NumberOfHinges(); i_hinge++ )
            {
                PointType hinge_coordiate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                printf("hinge %lu, block_Id %lu, (%lg, %lg, %lg)\n", 
                    i_hinge,
                    p_overset_condition->MeshBlockId(),
                    hinge_coordiate[0],
                    hinge_coordiate[1],
                    hinge_coordiate[2] );

                //
                HingeKey hinge_key{condition_id,i_hinge};
                HingeToAssignmentVectorMap::iterator it_assignment_key = hinge_to_assignments_map.find(hinge_key);

                if( it_assignment_key == hinge_to_assignments_map.end() )
                {
                    std::cout<<__func__<<"wrong! hinge_to_assignments_map"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                std::vector<PointSearchAssignmentKey> assignment_key_vector = it_assignment_key->second;

                int num_found = 0;
                for( const auto & assignment_key : assignment_key_vector )
                {
                    //
                    auto it_donor_info = donor_info_map.find(assignment_key);

                    if( it_donor_info == donor_info_map.end() )
                    {
                        std::cout<<__func__<<"wrong! donor_info_map"<<std::endl;
                        exit(EXIT_FAILURE);
                    }

                    //
                    PointSearchOutput donor_info = it_donor_info->second;
                    if( donor_info.mFound )
                    {
                        //dangerous write access to hinge
                        HingeDonorInfo & r_hinge_donor_info =  p_overset_condition->rHingeDonorInfo(i_hinge);

                        r_hinge_donor_info.mModelPartId              = donor_info.mModelPartId;
                        r_hinge_donor_info.mElementId                = donor_info.mElementId;
                        r_hinge_donor_info.mNodesId                  = donor_info.mNodesId;
                        r_hinge_donor_info.mBarycentricCoordinate[0] = donor_info.mBarycentricCoordinate[0];
                        r_hinge_donor_info.mBarycentricCoordinate[1] = donor_info.mBarycentricCoordinate[1];
                        r_hinge_donor_info.mBarycentricCoordinate[2] = donor_info.mBarycentricCoordinate[2];

                        printf("donor %lu (%lg, %lg, %lg), (%lg, %lg, %lg), found %d, distance %.10e \n", 
                            donor_info.mMeshBlockId,
                            donor_info.mInterpolatedCoordinate[0],
                            donor_info.mInterpolatedCoordinate[1],
                            donor_info.mInterpolatedCoordinate[2],
                            donor_info.mBarycentricCoordinate[0],
                            donor_info.mBarycentricCoordinate[1],
                            donor_info.mBarycentricCoordinate[2],
                            donor_info.mFound,
                            donor_info.mDistance );

                            num_found++;
                    }
                }

                printf("hinge found %d\n",num_found);

                if( num_found <= 0 )
                {
                    std::cout<<__func__<<"wrong! num_found"<<std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    void GetOversetConditionsDonorEquationsId()
    {
        InterpolateHingesDonorData();
    }

    void InterpolateHingesDonorData()
    {
        //
        struct HingeKey
        {
            OversetCondition::IndexType mConditionId;
            std::size_t mHingLocalId;

            bool operator< (const HingeKey & r_other ) const
            {
                if( mConditionId < r_other.mConditionId )
                    return true;
                else if ( mConditionId > r_other.mConditionId )
                    return false;
                else if ( mHingLocalId < r_other.mHingLocalId )
                    return true;

                return false;
            }
        };

        using HingeToAssignmentVectorMap = std::map<HingeKey,std::vector<InterpolationAssignmentKey>>;

        HingeToAssignmentVectorMap hinge_to_assignments_map;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            OversetCondition::IndexType condition_id = p_overset_condition->GetId();

            //add interpolation assignment
            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->NumberOfHinges(); i_hinge++ )
            {
                InterpolationInput interpolation_input;

                //dangerous write access to hinge
                HingeDonorInfo & r_hinge_donor_info =  p_overset_condition->rHingeDonorInfo(i_hinge);

                const std::size_t donor_model_part_id = r_hinge_donor_info.mModelPartId;

                interpolation_input.mElementId                = r_hinge_donor_info.mElementId;
                interpolation_input.mNodesId                  = r_hinge_donor_info.mNodesId;
                interpolation_input.mBarycentricCoordinate[0] = r_hinge_donor_info.mBarycentricCoordinate[0];
                interpolation_input.mBarycentricCoordinate[1] = r_hinge_donor_info.mBarycentricCoordinate[1];
                interpolation_input.mBarycentricCoordinate[2] = r_hinge_donor_info.mBarycentricCoordinate[2];
                
                //very bad!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //  assume interpolator has the same contractor key as dummy_model_part_holder
                InterpolatorKey r_interpolator_key = mModelPartIdToKey[donor_model_part_id];

                //add interpolation assignment
                InterpolationAssignmentKey interpolation_assignment_key = mpInterpolationMethod->AddInterpolationAssignment( r_interpolator_key, interpolation_input );

                //associate search assigment with hinge
                HingeKey hinge_key{condition_id,i_hinge};
                hinge_to_assignments_map[hinge_key].push_back(interpolation_assignment_key);
            }
        }

        //execute search
        mpInterpolationMethod->ExecuteAllInterpolationAssignments();

        //get search result
        std::vector<InterpolationMethod::InterpolationAssignmentOutputData> interpolation_output_data_vector;
        mpInterpolationMethod->GetInterpolationResults( interpolation_output_data_vector );


        //get search result mapped by assignment key
        std::map<InterpolationAssignmentKey,InterpolationOutput,InterpolationAssignmentKey::LessThanComparator> interpolation_output_map;

        for( const auto & r_interpolation_output_data : interpolation_output_data_vector )
        {
            interpolation_output_map[r_interpolation_output_data.GetAssignmentKey()] = r_interpolation_output_data.GetData();
        }

        // get result for hinges_data
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetConditions.ptr_begin(); it_p_condition != mOversetConditions.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition * p_overset_condition = dynamic_cast<OversetCondition *> ((* it_p_condition).get());

            OversetCondition::IndexType condition_id = p_overset_condition->GetId();

            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->NumberOfHinges(); i_hinge++ )
            {
                PointType hinge_coordiate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                printf("hinge %lu, block_Id %lu, (%lg, %lg, %lg)\n", 
                    i_hinge,
                    p_overset_condition->MeshBlockId(),
                    hinge_coordiate[0],
                    hinge_coordiate[1],
                    hinge_coordiate[2] );

                //
                HingeKey hinge_key{condition_id,i_hinge};
                HingeToAssignmentVectorMap::iterator it_assignment_key = hinge_to_assignments_map.find(hinge_key);

                if( it_assignment_key == hinge_to_assignments_map.end() )
                {
                    std::cout<<__func__<<"wrong! hinge_to_assignments_map"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                std::vector<InterpolationAssignmentKey> assignment_key_vector = it_assignment_key->second;

                for( const auto & r_assignment_key : assignment_key_vector )
                {
                    //
                    auto it_interpolation_output = interpolation_output_map.find(r_assignment_key);

                    if( it_interpolation_output == interpolation_output_map.end() )
                    {
                        std::cout<<__func__<<"wrong! hinge_donor_data_map"<<std::endl;
                        exit(EXIT_FAILURE);
                    }

                    //
                    InterpolationOutput & r_interpolation_output = it_interpolation_output->second;

                    HingeDonorData & r_hinge_donor_data = p_overset_condition->rHingeDonorData(i_hinge);

                    r_hinge_donor_data.mEquationsId = r_interpolation_output.mEquationsId;
                    r_hinge_donor_data.mNs          = r_interpolation_output.mNs;
                    r_hinge_donor_data.mDNsDXs      = r_interpolation_output.mDNsDXs;
                    r_hinge_donor_data.mTemperature = r_interpolation_output.mTemperature;
                    r_hinge_donor_data.mCoordinate  = r_interpolation_output.mCoordinate;

                    printf("donor (%lg, %lg, %lg), temp %lg\n", 
                        r_hinge_donor_data.mCoordinate[0],
                        r_hinge_donor_data.mCoordinate[1],
                        r_hinge_donor_data.mCoordinate[2],
                        r_hinge_donor_data.mTemperature );
                }
            }
        }
    }
    
private:
    const ModelPart & mrModelPart;
    OversetCommunicator mOversetCommunicator;
    DummyModelPartHolder mDummyModelPartHolder;
    DummyModelPartHolderManager mDummyModelPartHolderManager;
    std::map<Key,std::size_t,Key::LessThanComparator> mModelPartKeyToId;
    std::map<std::size_t,Key> mModelPartIdToKey;
    PointSearchMethod * mpPointSearchMethod;
    InterpolationMethod * mpInterpolationMethod;
    ModelPart::ConditionsContainerType mOversetConditions;
};

}//namespace OverserAssembly
}//namespace Kratos
#endif