#pragma once
#include"DistributedKey.h"

namespace DistributedAssignment
{

template<typename TLocationType>
class DistributedKeyIssuer
{
public:
    typedef DistributedKey<TLocationType> Key;

    DistributedKeyIssuer()
    :   mNumKeyIssued{0}
    {}

    ~DistributedKeyIssuer()
    {}

    void Clear()
    {
        mNumKeyIssued = 0;
    }

    Key IssueNewKey()
    {
        Key key = { TLocationType::Here(), mNumKeyIssued };
        mNumKeyIssued++;
        return key;
    }

private:
    std::size_t mNumKeyIssued;
};

}