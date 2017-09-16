#pragma once
#include <stdio.h>
#include <string>
#include "Point.h"
#include "Donor.h"


namespace SplitTreeSearch
{
#include "steForCppUser.hpp"
}

namespace DistributedPointSearcher
{

template<typename TContractorKeyType>
class SteSearcher
{
public:
    using Point = DistributedPointSearcher::Point ;
    using Donor = DistributedPointSearcher::Donor ;

    SteSearcher() = delete;

    SteSearcher( const std::vector<double> & r_nodes_coordinate, const std::vector<int> & r_elements_to_nodes, const int num_node, const int num_element )
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

        mSteHandle = SplitTreeSearch::steNew( mpCrd, mpCnn, PRM_TPL_4TET, num_element, SplitTreeSearch::TRUE )}
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
    
    void Execute( const Point & r_point, Donor & r_donor )
    {
        SplitTreeSearch::Real coordinate[3] = { r_point.mCoordinate[0], r_point.mCoordinate[1], r_point.mCoordinate[2] };
        SplitTreeSearch::Real local_coordinate[3];
        SplitTreeSearch::Real distance = 1e100;
    
        int element_id = SplitTreeSearch::steFindElemNextHeap( mSteHandle, coordinate , local_coordinate, & distance );
    
        r_donor.mElementId = element_id;
        r_donor.mLocalCoordinate[0] = (double) local_coordinate[0];
        r_donor.mLocalCoordinate[1] = (double) local_coordinate[1];
        r_donor.mLocalCoordinate[2] = (double) local_coordinate[2];
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

    std::string mName;
    TContractorKeyType mKey;

    SplitTreeSearch::Real * const mpCrd;
    SplitTreeSearch::Integer * const mpCnn;
    const SplitTreeSearch::Integer mNumNode;
    const SplitTreeSearch::Integer mNumElement;
    SplitTreeSearch::SteHd const mSteHandle;
    
    friend class DataUtility::DataPrinter;
};

}