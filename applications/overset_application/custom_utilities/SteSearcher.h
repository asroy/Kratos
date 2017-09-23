#if !defined(KRATOS_OVERSET_STESEARCHER_H_INCLUDED )
#define  KRATOS_OVERSET_STESEARCHER_H_INCLUDED

#include <stdio.h>
#include <string>
#include <vector>

#include "overset_application.h"

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
        const std::size_t mesh_block_id,
        const std::size_t num_node,
        const std::size_t num_element
    )
        :   mName(),
            mKey(),
            mpCrd{p_crd},
            mpCnn{p_cnn},
            mNodesId{local_id_to_node_id},
            mMeshBlockId{mesh_block_id},
            mNumNode{num_node},
            mNumElement{num_element},
            mSteHandle{SplitTreeSearch::steNew( mpCrd, mpCnn, PRM_TPL_4TET, mNumElement, SplitTreeSearch::TRUE )}
    {}

    virtual ~SteSearcher()
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
        using GlobalToLocal = std::map<std::size_t,std::size_t>;

        class BlockData
        {
            public:
            BlockData( const std::size_t num_node = 0, const std::size_t num_element = 0 )
                :   mNumberOfNode(num_node),
                    mNumberOfElement(num_element)
            {}

            std::size_t mNumberOfNode;
            std::size_t mNumberOfElement;
            std::vector<std::size_t> mCnn;
            std::vector<double> mCrd;
            GlobalToLocal mGlobalToLocal;
            std::vector<std::size_t> mLocalToGlobal;
        };

        using BlockDataMap = std::map<std::size_t, BlockData> ;

        //
        BlockDataMap blocks;

        const ModelPart::ElementsContainerType & r_elements_pointer = const_cast<ModelPart &>(r_model_part).Elements();
        
        //loop over element
        for( ModelPart::ElementsContainerType::ptr_const_iterator it_p_element = r_elements_pointer.ptr_begin(); it_p_element != r_elements_pointer.ptr_end(); it_p_element = std::next(it_p_element) )
        {
            const std::size_t block_id = (* it_p_element)->GetValue(BLOCK_ID);

            //if this block exist or not
            BlockDataMap::iterator it_block = blocks.find(block_id);

            //add a new block
            if( it_block == blocks.end() )
                it_block = blocks.insert( blocks.begin(), BlockDataMap::value_type {block_id, BlockData()} );

            BlockData & r_block = it_block->second;

            //add element
            r_block.mNumberOfElement = r_block.mNumberOfElement + 1;

            //loop over nodes
            Element::GeometryType r_nodes = (* it_p_element)->GetGeometry();

            for( std::size_t i = 0; i < r_nodes.size(); i++ )
            {
                const Element::NodeType & r_node = r_nodes[i];

                std::size_t global_id = r_node.GetId();

                GlobalToLocal::iterator it_local_id = r_block.mGlobalToLocal.find(global_id);

                std::size_t local_id;

                //add a new node
                if( it_local_id == r_block.mGlobalToLocal.end() )
                {
                    local_id = r_block.mNumberOfNode;

                    r_block.mNumberOfNode = r_block.mNumberOfNode + 1;
                    r_block.mLocalToGlobal.push_back(global_id);
                    r_block.mGlobalToLocal.insert( r_block.mGlobalToLocal.begin(), GlobalToLocal::value_type {global_id, local_id} );

                    r_block.mCrd.push_back( r_node.X() );
                    r_block.mCrd.push_back( r_node.Y() );
                    r_block.mCrd.push_back( r_node.Z() );
                }
                else
                {
                    local_id = it_local_id->second;
                }

                //add cnn
                r_block.mCnn.push_back(local_id);
            }
        }

        //create SteSearcher(s)
        r_point_searchers_pointer.clear();

        for ( BlockDataMap::const_iterator it_block = blocks.begin(); it_block != blocks.end(); it_block = std::next(it_block) )
        {
            const std::size_t block_id = it_block->first;
            const BlockData & r_block = it_block->second;

            const std::size_t num_node = r_block.mNumberOfNode;
            const std::size_t num_element = r_block.mNumberOfElement;

            double * const p_crd = new double [3*num_node];
            int * const p_cnn = new int [4*num_element];

            for( std::size_t i = 0; i < num_node; i++ )
            {
                p_crd[3*i]   = r_block.mCrd[3*i];
                p_crd[3*i+1] = r_block.mCrd[3*i+1];
                p_crd[3*i+2] = r_block.mCrd[3*i+2];
            }

            for( std::size_t i = 0; i < num_element; i++ )
            {
                p_cnn[4*i]   = r_block.mCnn[4*i];
                p_cnn[4*i+1] = r_block.mCnn[4*i+1];
                p_cnn[4*i+2] = r_block.mCnn[4*i+2];
                p_cnn[4*i+3] = r_block.mCnn[4*i+3];
            }

            r_point_searchers_pointer.push_back(new SteSearcher{ p_crd, p_cnn, r_block.mLocalToGlobal, block_id, num_node, num_element });

            std::cout<<__func__<<"block_id: " <<block_id<<", num_node: "<<num_node<<", num_element: "<<num_element<<std::endl;
            std::cout<<__func__<<"block_id: " <<block_id<<", size mCrd: "<<r_block.mCrd.size()<<", size mCnn: "<<r_block.mCnn.size()<<", size mLocalToGlobal: "<<r_block.mLocalToGlobal.size()<<", size mGlobalToLocal: "<<r_block.mGlobalToLocal.size()<<std::endl;
        }
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
    const std::size_t mMeshBlockId;
    const SplitTreeSearch::Integer mNumNode;
    const SplitTreeSearch::Integer mNumElement;

    const SplitTreeSearch::SteHd mSteHandle;
    
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
