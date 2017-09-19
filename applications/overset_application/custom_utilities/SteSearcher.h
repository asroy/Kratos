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
        const std::vector<double> & r_nodes_coordinate,
        const std::vector<int> & r_elements_to_nodes,
        const std::vector<std::size_t> & r_nodes_equation_id,
        const std::size_t num_node,
        const std::size_t num_element
    )
        :   mName(),
            mKey(),
            mpCrd{new double [r_nodes_coordinate.size()]},
            mpCnn{new int [r_elements_to_nodes.size()]},
            mNumNode{num_node},
            mNumElement{num_element},
            mSteHandle{nullptr}
    {
        for(std::size_t i = 0; i < r_nodes_coordinate.size(); i++ )
            mpCrd[i] = r_nodes_coordinate[i];

        for(std::size_t i = 0; i < r_elements_to_nodes.size(); i++ )
            mpCnn[i] = r_elements_to_nodes[i];

        mSteHandle = SplitTreeSearch::steNew( mpCrd, mpCnn, PRM_TPL_4TET, mNumElement, SplitTreeSearch::TRUE );

        mNodesEquationId.reserve(num_node);
        for( std::size_t i = 0; i < num_node; i++ )
            mNodesEquationId.push_back(r_nodes_equation_id[i]);
    }

    ~SteSearcher()
    {
        SplitTreeSearch::steFree( mSteHandle );
    }

    void SetName( const std::string & r_name )
    { mName = r_name; }

    std::string Name() const
    { return mName; }
    
    void SetKey(const SteSearcherKey key)
    { mKey = key; }
    
    SteSearcherKey Key() const
    { return mKey; }

    static void BuildFromModelPart(const ModelPart & r_model_part, std::vector<SteSearcher *> & r_point_searchers_pointer)
    {
        r_point_searchers_pointer.clear();

        //generate cnn, crd
        auto r_nodes = r_model_part.GetMesh().Nodes();
        auto r_elements = r_model_part.GetMesh().Elements();

        std::size_t num_node = r_nodes.size();
        std::size_t num_element = r_elements.size();

        double * p_crd = new double [3*num_node];
        int * p_cnn = new int [4*num_element];
        std::vector<std::size_t> node_local_to_equation_id(num_node);
        std::map<std::size_t,int> node_equation_to_local_id;

        {
            std::size_t i = 0;
            for( typename NodesArrayType::ptr_iterator it = r_model_part.GetMesh().Nodes().begin(); it ! = r_model_part.GetMesh().Nodes().end(); ++it )
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
            std::size_t i = 0;
            for( auto & r_element : r_model_part.GetMesh().Elements() )
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
        r_point_searchers_pointer.push_back(new SteSearcher( p_crd, p_ccn, node_local_to_equation_id, num_node, num_element ));
    }

    void Execute( const Coordinate & r_coordinate, DonorInfo & r_donor_info )
    {
        SplitTreeSearch::Real coordinate[3] = { r_coordinate.mCoordinate[0], r_coordinate.mCoordinate[1], r_coordinate.mCoordinate[2] };
        SplitTreeSearch::Real barycentric_coordinate[3];
        SplitTreeSearch::Real distance = 1e100;
    
        int element_id = SplitTreeSearch::steFindElemNextHeap( mSteHandle, coordinate , barycentric_coordinate, & distance );
        
        r_donor_info.mDonorNodesEquationId.clear();
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+1]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+2]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+3]] );

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
    const SplitTreeSearch::Integer mNumNode;
    const SplitTreeSearch::Integer mNumElement;
    SplitTreeSearch::SteHd mSteHandle;
    
    std::vector<std::size_t> mNodesEquationId;
    
    friend class DistributedAssignment::DataUtility::DataPrinter;
};

}//namespace OversetAssembly
}//namespace Kratos
#endif
