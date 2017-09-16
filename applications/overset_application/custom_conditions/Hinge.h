
template<std::size_t TDimension,
         typename THingeDataType,
         typename TDataType = double,
         typename TWeightType = double>
class Hinge : public IntegrationPoint<TDimension, TDataType, TWeightType>
{
private:
    using EquationIdVector = std::vector<std::size_t>;

public:
    ~Hinge() override;

private:
  	// A private default constructor necessary for serialization  
    Hinge();

    std::size_t mNumDonorNode;
    EquationIdVector mDonorNodesEquationIds;
    std::vector<THingeDataType> mDonorsData;

friend class Serializer;
};