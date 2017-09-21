#pragma once
#include<vector>

namespace DataUtility
{

//profile basic type
#define DATAPROFILE_PROFILE_BASIC_TYPE(TDATATYPE)     \
void Profile( const TDATATYPE & r_data )              \
{                                                     \
    this->SetIsTrivial(false);                        \
}

class DataProfile
{
public:
    DataProfile()
    :   mIsTrivial{true},
        mIsFromSender{false},
        mBufferContentSize{0}
    {}

    ~DataProfile()
    {}

    //profile basic data type
    DATAPROFILE_PROFILE_BASIC_TYPE(bool)
    DATAPROFILE_PROFILE_BASIC_TYPE(char)
    DATAPROFILE_PROFILE_BASIC_TYPE(int)
    DATAPROFILE_PROFILE_BASIC_TYPE(long)
    DATAPROFILE_PROFILE_BASIC_TYPE(unsigned int)
    DATAPROFILE_PROFILE_BASIC_TYPE(unsigned long)
    DATAPROFILE_PROFILE_BASIC_TYPE(float)
    DATAPROFILE_PROFILE_BASIC_TYPE(double)

    //profile std::vector
    template<typename TDataType>
    void Profile ( const std::vector<TDataType> & r_vector )
    {
        if( r_vector.size() == 0 )
            this->SetIsTrivial(true);
        else
            this->SetIsTrivial(false);
    }

    //profile user defined data type
    template<typename TDataType>
    void Profile ( const TDataType & r_data )
    {
        r_data.Profile(*this);
    }

    static DataProfile Default()
    {
        return DataProfile{true, true, 0};
    }

    DataProfile & SetIsTrivial( const bool is_trivial )
    {
        mIsTrivial = is_trivial;
        return *this;
    }

    DataProfile & SetIsFromSender( const bool is_from_sender )
    {
        mIsFromSender = is_from_sender;
        return *this;
    }

    DataProfile & SetBufferContentSize(const std::size_t buffer_size)
    {
        mBufferContentSize = buffer_size;
        return *this;
    }

    bool GetIsTrivial() const
    { return mIsTrivial; }

    bool GetIsFromSender() const
    { return mIsFromSender; }

    std::size_t GetBufferContentSize() const
    {
        if( mBufferContentSize < sizeof(Default()) )
        {
            std::cout<<__func__<<": mBufferContentSize too small! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }
        return mBufferContentSize;
     }

private:
    DataProfile( const bool is_trivial, const bool is_from_sender, const std::size_t buffer_content_size )
    :   mIsTrivial{is_trivial},
        mIsFromSender{is_from_sender},
        mBufferContentSize{buffer_content_size}
    {}

    bool mIsTrivial;
    bool mIsFromSender;
    std::size_t mBufferContentSize;
};

}