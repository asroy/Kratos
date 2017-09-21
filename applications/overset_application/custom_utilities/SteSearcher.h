#if !defined(KRATOS_OVERSET_STESEARCHER_H_INCLUDED )
#define  KRATOS_OVERSET_STESEARCHER_H_INCLUDED

#include <stdio.h>
#include <string>
#include <vector>

#include "DistributedAssignment.h"
#include "Coordinate.h"
#include "DonorInfo.h"

namespace Kratos
{
namespace OversetAssembly
{

namespace SplitTreeSearch
{
#include "SplitTreeSearch/steForCppUser.hpp"
}//namespace SplitTreeSearch


class SteSearcher
{
public:
    using Location = DistributedAssignment::Communication::MpiLocation;
    using SteSearcherKey = DistributedAssignment::DistributedAssignment::DistributedKey<Location>;

public:
    SteSearcher() = delete;

    SteSearcher
    ( 
        double * const p_crd,
        int * const p_cnn,
        const std::vector<std::size_t> local_id_to_node_id,
        const std::size_t num_node,
        const std::size_t num_element
    )
        :   mName(),
            mKey(),
            mpCrd{p_crd},
            mpCnn{p_cnn},
            mNodesId{local_id_to_node_id},
            mNumNode{num_node},
            mNumElement{num_element},
            mSteHandle{SplitTreeSearch::steNew( mpCrd, mpCnn, PRM_TPL_4TET, mNumElement, SplitTreeSearch::TRUE )}
    {}

    ~SteSearcher()
    {
        SplitTreeSearch::steFree( mSteHandle );//mpCrd, mpCnn will be deleted
    }

    void SetName( const std::string & r_name )
    { mName = r_name; }

    std::string Name() const
    { return mName; }
    
    void SetKey(const SteSearcherKey key)
    { mKey = key; }
    
    SteSearcherKey Key() const
    { return mKey; }

    static void BuildPointSearchersFromModelPart(const ModelPart & r_model_part, std::vector<SteSearcher *> & r_point_searchers_pointer)
    {
        r_point_searchers_pointer.clear();

        //generate cnn, crd
        std::size_t num_node = r_model_part.NumberOfNodes();
        std::size_t num_element = r_model_part.NumberOfElements();

        double * p_crd = new double [3*num_node];
        int * p_cnn = new int [4*num_element];
        std::vector<std::size_t> local_id_to_node_id(num_node);
        std::map<std::size_t,std::size_t> node_id_to_local_id;

        {
            const ModelPart::NodesContainerType & r_nodes_pointer = const_cast<ModelPart &>(r_model_part).Nodes();

            std::size_t i = 0;

            for( ModelPart::NodesContainerType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
            {
                p_crd[3*i]   = (* it_p_node)->X();
                p_crd[3*i+1] = (* it_p_node)->Y();
                p_crd[3*i+2] = (* it_p_node)->Z();
                i++;

                std::size_t node_id = (* it_p_node)->GetId();

                local_id_to_node_id[i] = node_id;
                node_id_to_local_id[node_id] = i;
            }
        }

        {
            const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(r_model_part).Elements();
            
            std::size_t i = 0;

            for( ModelPart::ElementsContainerType::ptr_const_iterator it_p_element = r_elements_pointer.ptr_begin(); it_p_element != r_elements_pointer.ptr_end(); it_p_element = std::next(it_p_element) )
            {
                const Element::GeometryType & r_nodes_pointer = (* it_p_element)->GetGeometry();

                assert( r_nodes_pointer.size() == 4 );

                for( Element::GeometryType::ptr_const_iterator it_p_node = r_nodes_pointer.ptr_begin(); it_p_node != r_nodes_pointer.ptr_end(); it_p_node = std::next(it_p_node) )
                {
                    p_cnn[i] = node_id_to_local_id[(* it_p_node)->GetId()];
                    i++;
                }
            }
        }
        
        //create SteSearcher
        r_point_searchers_pointer.push_back(new SteSearcher( p_crd, p_cnn, local_id_to_node_id, num_node, num_element ));
    }

    void Execute( const Coordinate & r_coordinate, DonorInfo & r_donor_info )
    {
        SplitTreeSearch::Real coordinate[3] = { r_coordinate.mCoordinate[0], r_coordinate.mCoordinate[1], r_coordinate.mCoordinate[2] };
        SplitTreeSearch::Real barycentric_coordinate[3];
        SplitTreeSearch::Real distance = 1e100;
    
        int element_id = SplitTreeSearch::steFindElemNextHeap( mSteHandle, coordinate , barycentric_coordinate, & distance );
        
        r_donor_info.mDonorNodesId.clear();
        r_donor_info.mDonorNodesId.push_back( mNodesId[mpCnn[4*element_id]] );
        r_donor_info.mDonorNodesId.push_back( mNodesId[mpCnn[4*element_id+1]] );
        r_donor_info.mDonorNodesId.push_back( mNodesId[mpCnn[4*element_id+2]] );
        r_donor_info.mDonorNodesId.push_back( mNodesId[mpCnn[4*element_id+3]] );

        r_donor_info.mBarycentricCoordinate[0] = (double) barycentric_coordinate[0];
        r_donor_info.mBarycentricCoordinate[1] = (double) barycentric_coordinate[1];
        r_donor_info.mBarycentricCoordinate[2] = (double) barycentric_coordinate[2];
    }

private:
    void Print( const DistributedAssignment::DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{SteSearcher: ";
        std::cout << "{Name: "<< mName <<"},",
        std::cout << "{ContractorKey: ";
        r_printer.Print(mKey);
        std::cout << "},";
        std::cout << "},";
    }

private:
    std::string mName;
    SteSearcherKey mKey;

    SplitTreeSearch::Real * const mpCrd;
    SplitTreeSearch::Integer * const mpCnn;
    const std::vector<std::size_t> mNodesId;
    const SplitTreeSearch::Integer mNumNode;
    const SplitTreeSearch::Integer mNumElement;
    const SplitTreeSearch::SteHd mSteHandle;
    
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
