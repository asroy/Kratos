#pragma once

namespace DistributedAssignment
{

template<typename TLocationType>
class DistributedKeyIssuer;

template<typename TLocationType>
class DistributedKey
{
public:
    DistributedKey()
    :   mBirthPlace(),
        mLocalKey{0}
    {}

    ~DistributedKey()
    {}

    static DistributedKey NoKey()
    {
        DistributedKey no_key = { TLocationType::NoWhere(), 0 };
        return no_key;
    }

    bool operator == (const DistributedKey & other) const
    {
        if( mBirthPlace != other.mBirthPlace )
            return false;
        else if( mLocalKey != other.mLocalKey )
            return false;
        else
            return true;
    }

    bool operator != (const DistributedKey & other) const
    { return !( *this == other ); }

    struct LessThanComparator
    {
        bool operator() ( const DistributedKey & a, const DistributedKey & b ) const
        {
            typename TLocationType::LessThanComparator location_less_than;

            if( location_less_than(a.mBirthPlace, b.mBirthPlace) )
                return true;
            else if( location_less_than(b.mBirthPlace, a.mBirthPlace) )
                return false;
            else if( a.mLocalKey < b.mLocalKey )
                return true;
            else
                return false;
        }
    };

private:
    DistributedKey( const TLocationType birth_place, const std::size_t local_key )
    :   mBirthPlace{birth_place},
        mLocalKey{local_key}
    {}

    void Save( DataUtility::Serializer & r_serializer ) const
    {
        r_serializer.Save(mBirthPlace);
        r_serializer.Save(mLocalKey);
    }

    void Load( DataUtility::Serializer & r_serializer )
    {
        r_serializer.Load(mBirthPlace);
        r_serializer.Load(mLocalKey);
    }

    void Profile( DataUtility::DataProfile & r_profile ) const
    {
        r_profile.SetIsTrivial(false);
    }

    void Print( const DataUtility::DataPrinter & r_printer ) const
    {
        std::cout << "{DistributedKey: ";
        r_printer.Print(mBirthPlace);
        r_printer.Print(mLocalKey);
        std::cout << "},";
    }

    TLocationType mBirthPlace;
    std::size_t mLocalKey;

    friend class DistributedKeyIssuer<TLocationType>;
    friend class DataUtility::Serializer;
    friend class DataUtility::DataProfile;
    friend class DataUtility::DataPrinter;
};

}