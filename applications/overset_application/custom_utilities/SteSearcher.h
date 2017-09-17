#pragma once
#include <stdio.h>
#include <string>
#include "Coordinate.h"
#include "DonorInfo.h"


namespace SplitTreeSearch
{
#include "steForCppUser.hpp"
}

namespace OversetAssembly
{

template<typename TContractorKeyType>
class SteSearcher
{
public:
    SteSearcher() = delete;

    SteSearcher
    ( 
        const std::vector<double> & r_nodes_coordinate,
        const std::vector<int> & r_elements_to_nodes,
        const std::vector<int> & r_nodes_equation_id,
        const int num_node,
        const int num_element
    )
        :   mName(),
            mKey(),
            mpCrd{nullptr},
            mpCnn{nullptr},
            mNumNode{num_node},
            mNumElement{num_element},
            mSteHandle{nullptr}
    {
        mpCrd = new double [r_nodes_coordinate.size()];
        mpCnn = new int [r_elements_to_nodes.size()];

        for(int i = 0; i < r_nodes_coordinate.size(); i++ )
            mpCrd[i] = r_nodes_coordinate[i];

        for(int i = 0; i < r_elements_to_nodes.size(); i++ )
            mpCnn[i] = r_elements_to_nodes[i];

        mSteHandle = SplitTreeSearch::steNew( mpCrd, mpCnn, PRM_TPL_4TET, num_element, SplitTreeSearch::TRUE )};

        mNodesEquationId.reserve(num_node);
        for( int i = 0; i < num_node; i++ )
            mNodesEquationId[i] = r_nodes_equation_id[i];
    }

    ~SteSearcher()
    {
        SplitTreeSearch::steFree( mSteHandle );
    }

    void SetName( const std::string & r_name )
    { mName = r_name; }

    std::string Name() const
    { return mName; }
    
    void SetKey(const TContractorKeyType key)
    { mKey = key; }
    
    TContractorKeyType Key() const
    { return mKey; }
    
    void Execute( const Coordinate & r_coordinate, DonorInfo & r_donor_info )
    {
        SplitTreeSearch::Real coordinate[3] = { r_coordinate.mCoordinate[0], r_coordinate.mCoordinate[1], r_coordinate.mCoordinate[2] };
        SplitTreeSearch::Real barycentric_coordinate[3];
        SplitTreeSearch::Real distance = 1e100;
    
        int element_id = SplitTreeSearch::steFindElemNextHeap( mSteHandle, coordinate , barycentric_coordinate, & distance );
        
        r_donor_info.mNumDonorNode = 4;

        r_donor_info.mDonorNodesEquationId.clear();
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+1]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+2]] );
        r_donor_info.mDonorNodesEquationId.push_back( mNodesEquationId[mpCnn[4*element_id+3]] );

        r_donor_info.mBarycentricCoordinate[0] = (double) local_coordinate[0];
        r_donor_info.mBarycentricCoordinate[1] = (double) local_coordinate[1];
        r_donor_info.mBarycentricCoordinate[2] = (double) local_coordinate[2];
    }

private:
    void Print( const DataUtility::DataPrinter & r_printer ) const
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
    TContractorKeyType mKey;

    SplitTreeSearch::Real * const mpCrd;
    SplitTreeSearch::Integer * const mpCnn;
    const SplitTreeSearch::Integer mNumNode;
    const SplitTreeSearch::Integer mNumElement;
    SplitTreeSearch::SteHd const mSteHandle;
    
    std::vector<std::size_t> mNodesEquationId;
    
    friend class DataUtility::DataPrinter;
};

}