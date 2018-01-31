#pragma once
#include<iostream>
#include<vector>

namespace DataUtility
{

//print basic type
#define DATAPRINTER_PRINT_BASIC_TYPE(TDATATYPE)         \
void Print( const TDATATYPE & r_data ) const            \
{                                                       \
    std::cout<<r_data<<",";                             \
}

class DataPrinter
{
  public:
    DataPrinter()
    {}

    ~DataPrinter()
    {}

    //print basic data type
    DATAPRINTER_PRINT_BASIC_TYPE(bool)
    DATAPRINTER_PRINT_BASIC_TYPE(char)
    DATAPRINTER_PRINT_BASIC_TYPE(int)
    DATAPRINTER_PRINT_BASIC_TYPE(long)
    DATAPRINTER_PRINT_BASIC_TYPE(unsigned int)
    DATAPRINTER_PRINT_BASIC_TYPE(unsigned long)
    DATAPRINTER_PRINT_BASIC_TYPE(float)
    DATAPRINTER_PRINT_BASIC_TYPE(double)

    //print std::string
    void Print ( const std::string & r_vector ) const
    {
        std::cout << r_vector;
    }

    //print std::vector
    template<typename TDataType>
    void Print ( const std::vector<TDataType> & r_vector ) const
    {
      std::cout<<"{std::vector, size "<<r_vector.size()<<", elements: "<<std::endl;

      for( const TDataType & r_data : r_vector )
      {
        Print(r_data);
        std::cout<<std::endl;
      }

      std::cout<<"}, "<<std::endl;
    }

    //print std::map
    template<typename TKeyType, typename TDataType>
    void Print ( const std::map<TKeyType,TDataType> & r_map ) const
    {
      std::cout<<"{std::map, size "<<r_map.size()<<", elements: "<<std::endl;

      for( auto it = r_map.begin(); it != r_map.end(); it = std::next(it) )
      {
        Print(it->first);
        Print(it->second);
        std::cout<<std::endl;
      }

      std::cout<<"}, "<<std::endl;
    }

    //print user defined data type
    template<typename TDataType>
    void Print ( const TDataType & r_data ) const
    {
      r_data.Print(*this);
    }
};

}
