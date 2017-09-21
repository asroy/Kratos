#pragma once
#include<iostream>
#include<vector>
#include<map>
#include<mpi.h>
#include"Serializer.h"
#include"DataProfile.h"
#include"MpiLocation.h"

namespace Communication
{

class MpiCommunicator
{
public:
    using Location = MpiLocation;
    using Peer = int;

    template<typename TDataType>
    using MapByLocationType = std::map<Location, TDataType, typename Location::LessThanComparator>;

    template<typename TDataType>
    using MapByPeerType = std::map<Peer, TDataType>;

    template<typename TDataType>
    using PairByLocationType = std::pair<const Location, TDataType>;

    template<typename TDataType>
    using PairByPeerType = std::pair<const Peer, TDataType>;

    MpiCommunicator()
    :   mMpiComm{MPI_COMM_WORLD},
        mPeerToLocation(),
        mLocationToPeer(),
        mPeers(),
        mLocations(),
        mSendSerializers(),
        mRecvSerializers()
    {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, & size);

        for ( int i = 0; i < size; i++ )
        {
            Location location(i);
            Peer peer = i;
            mPeerToLocation[i] = location;
            mLocationToPeer[location] = peer;
            mPeers.push_back(peer);
            mLocations.push_back(location);
        }
    }

    MpiCommunicator(const MPI_Comm & mpi_comm)
    :   mMpiComm{mpi_comm},
        mPeerToLocation(),
        mLocationToPeer(),
        mPeers(),
        mLocations(),
        mSendSerializers(),
        mRecvSerializers()
    {
        int size;
        MPI_Comm_size(mMpiComm, & size);

        int * world_ranks = new int[size];

        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, & world_rank);

        MPI_Allgather( & world_rank, 1, MPI_INT, world_ranks, 1, MPI_INT, mMpiComm );

        for ( int i = 0; i < size; i++ )
        {
            Location location(world_ranks[i]);
            Peer peer = i;
            mPeerToLocation[i] = location;
            mLocationToPeer[location] = peer;
            mPeers.push_back(peer);
            mLocations.push_back(location);
        }

        delete[] world_ranks;
    }

    ~MpiCommunicator()
    {}

    Peer Myself() const
    {
        int rank;
        MPI_Comm_rank(mMpiComm, & rank);
        return rank;
    }

    const std::vector<Peer> & Peers() const
    { return mPeers; }

    const std::vector<Location> & Locations() const
    { return mLocations; }

    Location PeerToLocation(const Peer & peer) const
    {
        const typename MapByPeerType<Location>::const_iterator it = mPeerToLocation.find(peer);

        if( it == mPeerToLocation.end() )
        {
            std::cout<<__func__<<"peer not exist! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }

    Peer LocationToPeer(const Location & location) const
    {
        const typename MapByLocationType<Peer>::const_iterator it = mLocationToPeer.find(location);

        if( it == mLocationToPeer.end() )
        {
            std::cout<<__func__<<"location not exist! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return it->second;
    }

    template<typename TDataType>
    void AllSendAllRecv( const MapByLocationType<TDataType> & r_send_datas, MapByLocationType<TDataType> & r_recv_datas, const int mpi_tag )
    {
        typedef PairByLocationType<TDataType> DataPair;
        typedef PairByLocationType<DataUtility::Serializer> SerializerPair;

        //clear
        mSendSerializers.clear();
        mRecvSerializers.clear();
        r_recv_datas.clear();

        //serializer r_send_datas
        for( const DataPair & r_send_data_pair: r_send_datas )
        {
            const Location send_to_location = r_send_data_pair.first;
            const TDataType & r_send_data = r_send_data_pair.second;

            DataUtility::Serializer & r_send_serializer = mSendSerializers[send_to_location];

            //save to serializers
            r_send_serializer.WriteBufferHeader(DataUtility::DataProfile::Default());
            std::size_t send_size = r_send_serializer.FreshSave(r_send_data);

            //update header
            DataUtility::DataProfile send_data_profile;
            send_data_profile.Profile(r_send_data);
            send_data_profile.SetIsFromSender(true).SetBufferContentSize(send_size);
            r_send_serializer.WriteBufferHeader(send_data_profile);
        }

        //send recv buffer
        AllSendAllRecvSerializer(mSendSerializers, mRecvSerializers, mpi_tag);

        //de-serializer r_recv_datas
        for( SerializerPair & r_recv_serializer_pair : mRecvSerializers )
        {
            Location recv_location = r_recv_serializer_pair.first;
            DataUtility::Serializer & r_recv_serializer = r_recv_serializer_pair.second;

            DataUtility::DataProfile recv_data_profile;
            r_recv_serializer.ReadBufferHeader(recv_data_profile);

            //check
            if( ! recv_data_profile.GetIsFromSender() )
            {
                std::cout<<__func__<<": recv serializer exit, but not from sender! exit"<<std::endl;
                exit(EXIT_FAILURE);
            }

            //check
            if( recv_data_profile.GetIsTrivial() )
            {
                std::cout<<__func__<<": recv serializer exit, but buffer header is marked trivial! exit"<<std::endl;
                exit(EXIT_FAILURE);
            }

            if( recv_data_profile.GetIsFromSender() )
            {
                TDataType & r_recv_data = r_recv_datas[recv_location];//this will insert a entry into r_recv_datas map
                r_recv_serializer.FreshLoad(r_recv_data);
            }
        }

        //clean up serializers
        mSendSerializers.clear();
        mRecvSerializers.clear();
    }

    template<typename TDataType>
    void AllGather( const TDataType & r_send_data, MapByLocationType<TDataType> & r_recv_datas, int mpi_tag )
    {
        typedef PairByLocationType<DataUtility::Serializer> SerializerPair;

        //clear r_recv_datas
        r_recv_datas.clear();

        //root location
        Location my_location = PeerToLocation(Myself());

        //serialize r_send_data
        {
            DataUtility::Serializer & r_send_serializer = mSendSerializers[my_location];

            //save to serializer
            r_send_serializer.WriteBufferHeader(DataUtility::DataProfile::Default());
            std::size_t send_size = r_send_serializer.FreshSave(r_send_data);

            //update header
            DataUtility::DataProfile send_data_profile;
            send_data_profile.Profile(r_send_data);
            send_data_profile.SetIsFromSender(true).SetBufferContentSize(send_size);
            r_send_serializer.WriteBufferHeader(send_data_profile);
        }

        //send recv buffer, recv ser
        AllGatherSerializer( mSendSerializers[my_location], mRecvSerializers, mpi_tag );

        //de-serializer r_recv_datas
        for( SerializerPair & r_recv_serializer_pair: mRecvSerializers )
        {
            Location recv_location = r_recv_serializer_pair.first;
            DataUtility::Serializer & r_recv_serializer = r_recv_serializer_pair.second;

            DataUtility::DataProfile recv_data_profile;
            r_recv_serializer.ReadBufferHeader(recv_data_profile);

            //check
            if( ! recv_data_profile.GetIsFromSender() )
            {
                std::cout<<__func__<<": recv serializer exit, but not from sender! exit"<<std::endl;
                exit(EXIT_FAILURE);
            }

            //check
            if( recv_data_profile.GetIsTrivial() )
            {
                std::cout<<__func__<<": recv serializer exit, but buffer header is marked trivial! exit"<<std::endl;
                exit(EXIT_FAILURE);
            }

            if( recv_data_profile.GetIsFromSender() )
            {
                TDataType & r_recv_data = r_recv_datas[recv_location];//this will insert a entry into r_recv_datas map
                r_recv_serializer.FreshLoad(r_recv_data);
            }
        }

        //clean up serializers
        mSendSerializers.clear();
        mRecvSerializers.clear();
    }

private:

    int PeerToRank(const Peer & peer) const
    {
        int size;
        MPI_Comm_size(mMpiComm, & size);

        int rank = peer;
        if( rank < 0 || rank >= size )
        {
            std::cout<<__func__<<"peer not exist! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return rank;
     }

    Peer RankToPeer(const int rank) const
    {
        int size;
        MPI_Comm_size(mMpiComm, & size);

        if( rank < 0 || rank >= size )
        {
            std::cout<<__func__<<"rank not exist! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        return rank;
    }

    void AllSendAllRecvSerializer( const MapByLocationType<DataUtility::Serializer> & r_send_serializer_map, MapByLocationType<DataUtility::Serializer> & r_recv_serializer_map, int mpi_tag )
    {
        typedef PairByLocationType<DataUtility::Serializer> SerializerPair;

        //clear recv serializers
        r_recv_serializer_map.clear();

        int mpi_size;
        MPI_Comm_size(mMpiComm, &mpi_size);

        int* send_sizes = new int[mpi_size];
        int* recv_sizes = new int[mpi_size];

        for( int i =0; i < mpi_size; i++ )
        {
            send_sizes[i] = 0;
            recv_sizes[i] = 0;
        }

        //message sizes
        for( const SerializerPair & r_send_serializer_pair : r_send_serializer_map )
        {
            Location send_to_location = r_send_serializer_pair.first;
            const DataUtility::Serializer & r_send_serializer = r_send_serializer_pair.second;

            DataUtility::DataProfile send_data_profile;
            r_send_serializer.ReadBufferHeader(send_data_profile);

            int i = PeerToRank(LocationToPeer(send_to_location));
            send_sizes[i] = send_data_profile.GetBufferContentSize();

            if ( send_data_profile.GetIsTrivial() )//don't want to send trivial data
                send_sizes[i] = 0;
        }

        //receive message size
        MPI_Alltoall( send_sizes, 1, MPI_INT, recv_sizes, 1, MPI_INT, mMpiComm );

        //reserve space and write dummy header in recv buffer
        for( int i = 0; i < mpi_size; i++ )
        {
            if(recv_sizes[i] > 0)
            {
                Location recv_from_location = PeerToLocation(RankToPeer(i));
                DataUtility::Serializer & r_recv_serializer = r_recv_serializer_map[recv_from_location];//this will insert new entry into map
                r_recv_serializer.IncreaseBufferSize(recv_sizes[i]);
                r_recv_serializer.WriteBufferHeader(DataUtility::DataProfile::Default().SetIsFromSender(false).SetIsTrivial(true));
            }
        }

        int num_event = 0;
        for ( int i = 0; i < mpi_size; i++ )
        {
            if( send_sizes[i] > 0 ) num_event++;
            if( recv_sizes[i] > 0 ) num_event++;
        }

        MPI_Request * reqs = new MPI_Request [num_event];
        MPI_Status * stats = new MPI_Status  [num_event];

        // send and recv buffer, trivial data will not be received
        num_event = 0;
        for( int i = 0; i < mpi_size; i++ )
        {
            if( recv_sizes[i] > 0 )
            {
                Location recv_from_location = PeerToLocation(RankToPeer(i));
                DataUtility::Serializer & r_recv_serializer = r_recv_serializer_map[recv_from_location];
                MPI_Irecv( r_recv_serializer.BufferPointer(), recv_sizes[i], MPI_CHAR, i, mpi_tag, mMpiComm, &(reqs[num_event]) );
                num_event++;
            }

            if( send_sizes[i] > 0 )
            {
                Location send_to_location = PeerToLocation(RankToPeer(i));
                const DataUtility::Serializer & r_send_serializer = mSendSerializers[send_to_location];
                MPI_Isend( r_send_serializer.BufferPointer(), send_sizes[i], MPI_CHAR, i, mpi_tag, mMpiComm, &(reqs[num_event]) );
                num_event++;
            }
        }

        // wait until all send recv are done
        MPI_Waitall( num_event, reqs, stats );

        delete [] send_sizes;
        delete [] recv_sizes;
        delete [] reqs;
        delete [] stats;
    }

    void AllGatherSerializer( const DataUtility::Serializer & r_send_serializer, MapByLocationType<DataUtility::Serializer>  & r_recv_serializer_map, int mpi_tag )
    {
        //clear recv serializers
        r_recv_serializer_map.clear();

        int mpi_size;
        MPI_Comm_size(mMpiComm, &mpi_size);

        int* send_sizes = new int[mpi_size];
        int* recv_sizes = new int[mpi_size];

        for( int i = 0; i < mpi_size; i++ )
        {
            send_sizes[i] = 0;
            recv_sizes[i] = 0;
        }

        //serialize send_data
        int send_size;
        DataUtility::DataProfile send_data_profile;
        r_send_serializer.ReadBufferHeader(send_data_profile);

        send_size = send_data_profile.GetBufferContentSize();

        if ( send_data_profile.GetIsTrivial() )//don't want to send trivial data
            send_size = 0;

        //
        for( int i = 0; i < mpi_size; i++ )
            send_sizes[i] = send_size;

        //receive message size
        MPI_Alltoall( send_sizes, 1, MPI_INT, recv_sizes, 1, MPI_INT, mMpiComm );

        //reserve space and write dummy header in recv buffer
        for( int i = 0; i < mpi_size; i++ )
        {
            if(recv_sizes[i] > 0)
            {
                Location recv_from_location = PeerToLocation(RankToPeer(i));
                DataUtility::Serializer & r_recv_serializer = r_recv_serializer_map[recv_from_location];

                r_recv_serializer.IncreaseBufferSize(recv_sizes[i]);
                r_recv_serializer.WriteBufferHeader( DataUtility::DataProfile::Default().SetIsFromSender(false).SetIsTrivial(true) );
            }
        }

        int num_event = 0;
        for ( int i = 0; i < mpi_size; i++ )
        {
            if( send_sizes[i] > 0 ) num_event++;
            if( recv_sizes[i] > 0 ) num_event++;
        }

        MPI_Request * reqs = new MPI_Request [num_event];
        MPI_Status * stats = new MPI_Status  [num_event];

        // send and recv buffer, trivial data will not be sent
        num_event = 0;
        for( int i = 0; i < mpi_size; i++ )
        {
            if( recv_sizes[i] > 0 )
            {
                Location recv_from_location = PeerToLocation(RankToPeer(i));
                DataUtility::Serializer & r_recv_serializer = r_recv_serializer_map[recv_from_location];

                MPI_Irecv( r_recv_serializer.BufferPointer(), recv_sizes[i], MPI_CHAR, i, mpi_tag, mMpiComm, &(reqs[num_event]) );
                num_event++;
            }

            if( send_size > 0 )
            {
                MPI_Isend( r_send_serializer.BufferPointer(), send_size, MPI_CHAR, i, mpi_tag, mMpiComm, &(reqs[num_event]) );
                num_event++;
            }
        }

        // wait until all send recv are done
        MPI_Waitall( num_event, reqs, stats );

        delete [] send_sizes;
        delete [] recv_sizes;
        delete [] reqs;
        delete [] stats;
    }

    MPI_Comm mMpiComm;
    MapByPeerType<Location> mPeerToLocation;
    MapByLocationType<Peer> mLocationToPeer;

    std::vector<Peer> mPeers;
    std::vector<Location> mLocations;

    MapByLocationType<DataUtility::Serializer> mSendSerializers;
    MapByLocationType<DataUtility::Serializer> mRecvSerializers;
};

}