
template<std::size_t TDimension,
         typename THingeData,
         typename TData = double,
         typename TWeight = double>
class Hinge : public IntegrationPoint<TDimension, TData, TWeight>
{
//type
private:
    using EquationIdVector = std::vector<std::size_t>;

//method
public:
    ~Hinge() override;

private:
  	// A private default constructor necessary for serialization  
    Hinge();

//member
private:
    EquationIdVector mDonorNodesEquationIds;
    std::vector<THingeData> mDonorsData;

friend class Serializer;
};