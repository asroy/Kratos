#pragma once

template<std::size_t TDimension,
         typename THingeData,
         typename TData = double,
         typename TWeight = double>
class Hinge : public IntegrationPoint<TDimension, TData, TWeight>
{
private:
    using Base = IntegrationPoint<TDimension, TData, TWeight>;
    
public:
    Hinge( const IntegrationPoint & r_integration_point )
        :   Base{r_integration_point}
    {}

    ~Hinge() override
    {}

private:
  	// A private default constructor necessary for serialization  
    Hinge()
        :   Base()
    {}

//member
private:
    std::vector<std::size_t> mDonorNodesEquationIds;
    std::vector<THingeData> mDonorsData;

friend class Serializer;
};