#pragma once
#include<iostream>
#include<mpi.h>

namespace Communication
{

//implementation of Location type for MpiCommunicator
class MpiLocation
{
public:
    MpiLocation()
    :   mWorldRank{-1}
    {}

    ~MpiLocation()
    {}

    static MpiLocation Here()
    {
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, & world_rank);
        return MpiLocation{world_rank};
    }

    static MpiLocation NoWhere()
    {
        return MpiLocation{-1};
    }

    MpiLocation Next( const int n ) const
    {
        return MpiLocation{mWorldRank + n};
    }

    bool operator == (const MpiLocation & other) const
    {
        if( mWorldRank != other.mWorldRank )
            return false;
        else
            return true;
    }

    bool operator != (const MpiLocation & other) const
    { return !( *this == other ); }

    struct LessThanComparator
    {
        bool operator() ( const MpiLocation & a, const MpiLocation & b ) const
        {
            return ( a.mWorldRank < b.mWorldRank );
        }
    };

private:
    MpiLocation( const int world_rank )
    :   mWorldRank{world_rank}
    {}

    void Save( DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mWorldRank);
    }

    void Load( DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mWorldRank);
    }

    void Profile( DataUtility::DataProfile & r_profile ) const
    {
       r_profile.SetIsTrivial(false);
    }

    void Print( const DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{MpiLocation: ";
        r_printer.Print(mWorldRank);
        std::cout << "},";
    }

    int mWorldRank;

    friend class MpiCommunicator;
    friend class DataUtility::Serializer;
    friend class DataUtility::DataProfile;
    friend class DataUtility::DataPrinter;
};

}