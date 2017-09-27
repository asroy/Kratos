#if !defined(KRATOS_OVERSET_ASSEMBLY_H_INCLUDED )
#define  KRATOS_OVERSET_ASSEMBLY_H_INCLUDED

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_conditions/OversetCondition3D.h"
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

public:
    OversetAssembler() = delete;

    OversetAssembler(const ModelPart & r_model_part)
        :   mrModelPart{r_model_part},
            mOversetCommunicator(),
            mDummyModelPartHolder(),
            mDummyModelPartHolderManager{mOversetCommunicator},
            mpPointSearchMethod{nullptr},
            mpInterpolationMethod{nullptr},
            mOversetCondition3Ds()
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
            }
        }

        std::cout<<__func__<<": size face_to_element_map: "<<face_to_element_map.size()<<std::endl;


        //loop over overset conditions
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
        std::size_t num_overset_condition = 0;
        std::size_t num_hinge = 0;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            if ( ! p_overset_condition )
            {
                //throw error please
                std::cout<<__func__<<"wrong! not OversetCondition3D *"<<std::endl;
                exit(EXIT_FAILURE);
            }

            p_overset_condition->GenerateHinges();

            num_overset_condition++;
            num_hinge += p_overset_condition->Hinge3Ds().size();
        }

        std::cout<<__func__<<": num_overset_condition "<<num_overset_condition<<", num_hinge "<<num_hinge<<std::endl;
    }

    void GenerateHingeDonorRelation()
    {
        //
        struct HingeId
        {
            OversetCondition3D::IndexType mConditionId;
            std::size_t mHingLocalId;

            bool operator< (const HingeId & r_other ) const
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

        using HingeToAssignmentVectorMap = std::map<HingeId,std::vector<PointSearchAssignmentKey>>;

        HingeToAssignmentVectorMap hinge_to_assignments_map;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            OversetCondition3D::IndexType condition_id = p_overset_condition->GetId();

            //add search assignment
            std::size_t condition_block_id = p_overset_condition->MeshBlockId();

            for ( const PointSearchMethod::PointSearcherKeySetMapByMeshBlockId::value_type & r_pair : mpPointSearchMethod->GlobalSearchersKeyForBlock() )
            {
                std::size_t searcher_block_id = r_pair.first;

                if( condition_block_id != searcher_block_id )
                {
                    for( const PointSearcherKey & r_searcher_key : r_pair.second )
                    {
                        for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->Hinge3Ds().size(); i_hinge++ )
                        {  
                            const Vector hinge_coordinate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                            //add search assignment
                            PointSearchAssignmentKey search_assignment_key = mpPointSearchMethod->AddSearch( r_searcher_key, {hinge_coordinate[0], hinge_coordinate[1], hinge_coordinate[2]} );

                            //associate search assigment with hinge
                            HingeId hinge_key{condition_id,i_hinge};
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
        std::map<PointSearchAssignmentKey,DonorInfo,PointSearchAssignmentKey::LessThanComparator> donor_info_map;

        for( const auto & r_donor_info_data : donor_info_data_vector )
        {
            donor_info_map[r_donor_info_data.GetAssignmentKey()] = r_donor_info_data.GetData();
        }

        // get donor result for hinges
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            OversetCondition3D::IndexType condition_id = p_overset_condition->GetId();

            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->Hinge3Ds().size(); i_hinge++ )
            {
                Vector hinge_coordiate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                printf("hinge %lu, block_Id %lu, (%lg, %lg, %lg)\n", 
                    i_hinge,
                    p_overset_condition->MeshBlockId(),
                    hinge_coordiate[0],
                    hinge_coordiate[1],
                    hinge_coordiate[2] );

                //
                HingeId hinge_key{condition_id,i_hinge};
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
                    DonorInfo donor_info = it_donor_info->second;
                    if( donor_info.mFound )
                    {
                        //dangerous write access to hinge
                        Hinge3D & r_hinge =  p_overset_condition->rHinge3D(i_hinge);

                        r_hinge.rDonorModelPartId()              = donor_info.mDonorModelPartId;
                        r_hinge.rDonorElementId()                = donor_info.mDonorElementId;
                        r_hinge.rDonorNodesId()                  = donor_info.mDonorNodesId;
                        r_hinge.rDonorBarycentricCoordinate()[0] = donor_info.mDonorBarycentricCoordinate[0];
                        r_hinge.rDonorBarycentricCoordinate()[1] = donor_info.mDonorBarycentricCoordinate[1];
                        r_hinge.rDonorBarycentricCoordinate()[2] = donor_info.mDonorBarycentricCoordinate[2];

                        printf("donor %lu (%lg, %lg, %lg), (%lg, %lg, %lg), found %d, distance %.10e \n", 
                            donor_info.mDonorMeshBlockId,
                            donor_info.mInterpolatedCoordinate.mCoordinate[0],
                            donor_info.mInterpolatedCoordinate.mCoordinate[1],
                            donor_info.mInterpolatedCoordinate.mCoordinate[2],
                            donor_info.mDonorBarycentricCoordinate[0],
                            donor_info.mDonorBarycentricCoordinate[1],
                            donor_info.mDonorBarycentricCoordinate[2],
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

    void InterpolateHingeData()
    {
        //
        struct HingeId
        {
            OversetCondition3D::IndexType mConditionId;
            std::size_t mHingLocalId;

            bool operator< (const HingeId & r_other ) const
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

        using HingeToAssignmentVectorMap = std::map<HingeId,std::vector<InterpolationAssignmentKey>>;

        HingeToAssignmentVectorMap hinge_to_assignments_map;

        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            OversetCondition3D::IndexType condition_id = p_overset_condition->GetId();

            //add interpolation assignment
            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->Hinge3Ds().size(); i_hinge++ )
            {
                InterpolationInput interpolation_input;

                //dangerous write access to hinge
                Hinge3D & r_hinge =  p_overset_condition->rHinge3D(i_hinge);

                interpolation_input.mElementId                = r_hinge.rDonorElementId();
                interpolation_input.mNodesId                  = r_hinge.rDonorNodesId();
                interpolation_input.mBarycentricCoordinate[0] = r_hinge.rDonorBarycentricCoordinate()[0];
                interpolation_input.mBarycentricCoordinate[1] = r_hinge.rDonorBarycentricCoordinate()[1];
                interpolation_input.mBarycentricCoordinate[2] = r_hinge.rDonorBarycentricCoordinate()[2];

                const std::size_t donor_model_part_id = r_hinge.rDonorModelPartId();
                
                //very bad!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                //  assume interpolator has the same contractor key as dummy_model_part_holder
                InterpolatorKey r_interpolator_key = mModelPartIdToKey[donor_model_part_id];

                //add interpolation assignment
                InterpolationAssignmentKey interpolation_assignment_key = mpInterpolationMethod->AddInterpolationAssignment( r_interpolator_key, interpolation_input );

                //associate search assigment with hinge
                HingeId hinge_key{condition_id,i_hinge};
                hinge_to_assignments_map[hinge_key].push_back(interpolation_assignment_key);
            }
        }

        //execute search
        mpInterpolationMethod->ExecuteAllInterpolationAssignments();

        //get search result
        std::vector<InterpolationMethod::InterpolationAssignmentOutputData> hinge_data_data_vector;
        mpInterpolationMethod->GetInterpolationResults( hinge_data_data_vector );


        //get search result mapped by assignment key
        std::map<InterpolationAssignmentKey,HingeData,InterpolationAssignmentKey::LessThanComparator> hinge_data_map;

        for( const auto & r_hinge_data_data : hinge_data_data_vector )
        {
            hinge_data_map[r_hinge_data_data.GetAssignmentKey()] = r_hinge_data_data.GetData();
        }

        // get result for hinges_data
        for( ModelPart::ConditionsContainerType::ptr_const_iterator it_p_condition = mOversetCondition3Ds.ptr_begin(); it_p_condition != mOversetCondition3Ds.ptr_end(); it_p_condition = std::next(it_p_condition) )
        {
            OversetCondition3D * p_overset_condition = dynamic_cast<OversetCondition3D *> ((* it_p_condition).get());

            OversetCondition3D::IndexType condition_id = p_overset_condition->GetId();

            for( std::size_t i_hinge = 0; i_hinge < p_overset_condition->Hinge3Ds().size(); i_hinge++ )
            {
                Vector hinge_coordiate = p_overset_condition->HingeGlobalCoordinate(i_hinge);

                printf("hinge %lu, block_Id %lu, (%lg, %lg, %lg)\n", 
                    i_hinge,
                    p_overset_condition->MeshBlockId(),
                    hinge_coordiate[0],
                    hinge_coordiate[1],
                    hinge_coordiate[2] );

                //
                HingeId hinge_key{condition_id,i_hinge};
                HingeToAssignmentVectorMap::iterator it_assignment_key = hinge_to_assignments_map.find(hinge_key);

                if( it_assignment_key == hinge_to_assignments_map.end() )
                {
                    std::cout<<__func__<<"wrong! hinge_to_assignments_map"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                std::vector<InterpolationAssignmentKey> assignment_key_vector = it_assignment_key->second;

                for( const auto & assignment_key : assignment_key_vector )
                {
                    //
                    auto it_hinge_data = hinge_data_map.find(assignment_key);

                    if( it_hinge_data == hinge_data_map.end() )
                    {
                        std::cout<<__func__<<"wrong! hinge_data_map"<<std::endl;
                        exit(EXIT_FAILURE);
                    }

                    //
                    HingeData hinge_data = it_hinge_data->second;

                    printf("donor (%lg, %lg, %lg) \n", 
                        hinge_data.mCoordinate[0],
                        hinge_data.mCoordinate[1],
                        hinge_data.mCoordinate[2] );
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
    ModelPart::ConditionsContainerType mOversetCondition3Ds;
};

}//namespace OverserAssembly
}//namespace Kratos
#endif