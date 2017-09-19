#pragma once
#include<iostream>
#include<vector>
#include<map>
#include"DistributedContractorManager.h"
#include"AssignmentData.h"

namespace DistributedAssignment{

template<   typename TAssignorType,
            typename TAssigneeType,
            typename TInputType,
            typename TOutputType,
            typename TCommunicatorType,
            template <typename TDummyType> class TDistributedContractorKeyIssuerType,
            template <typename TDummyType> class TDistributedAssignmentKeyIssuerType >
class DistributedAssignmentManager
{
public:
    using Location = typename TCommunicatorType::Location ;
    using ContractorKeyIssuer = TDistributedContractorKeyIssuerType<Location> ;
    using ContractorKey = typename ContractorKeyIssuer::Key ;
    using AssignmentKeyIssuer = TDistributedAssignmentKeyIssuerType<Location> ;
    using AssignmentKey = typename AssignmentKeyIssuer::Key ;

    template<typename TContractorType>
    using ContractorManagerType = DistributedContractorManager<TContractorType, TCommunicatorType, TDistributedContractorKeyIssuerType>;

    template<typename TDataType>
    using AssignmentDataType = AssignmentData<ContractorKey,AssignmentKey,TDataType>;

    template<typename TDataType>
    using AssignmentDataVectorType = std::vector<AssignmentDataType<TDataType>>;

    template<typename TKeyType, typename TDataType>
    using AssignmentDataVectorMapType = std::map<TKeyType, AssignmentDataVectorType<TDataType>, typename TKeyType::LessThanComparator>;

    template<typename TKeyType, typename TDataType>
    using AssignmentDataVectorPairType = std::pair<const TKeyType,  AssignmentDataVectorType<TDataType>>;

    DistributedAssignmentManager() = delete;

    DistributedAssignmentManager
    (   TCommunicatorType & r_communicator,
        ContractorManagerType<TAssignorType> & r_assignor_manager,
        ContractorManagerType<TAssigneeType> & r_assignee_manager
    )
        :   mpCommunicator {& r_communicator},
            mpAssignorManager {& r_assignor_manager},
            mpAssigneeManager {& r_assignee_manager},
            mAssignmentKeyIssuer ()
    {}

    virtual ~DistributedAssignmentManager()
    {
        mpCommunicator = nullptr;
        mpAssignorManager = nullptr;
        mpAssigneeManager = nullptr;
    }

    void ClearAllAssignments()
    {
        mAssignmentKeyIssuer.Clear();

        mAssignorInputs.clear();
        mAssigneeInputs.clear();
        mWorkUnitInputs.clear();
        mAssignorOutputs.clear();
        mAssigneeOutputs.clear();
        mWorkUnitOutputs.clear();
    }

    AssignmentKey AddAssignment( const ContractorKey assignor_key, const ContractorKey assignee_key, const TInputType input_data )
    {
        typedef AssignmentDataType<TInputType> AssignmentInput;

        AssignmentKey assignment_key = mAssignmentKeyIssuer.IssueNewKey();

        AssignmentInput assignment_input = { assignor_key, assignee_key, assignment_key, input_data };

        Location assignee_location = mpAssigneeManager->FindGlobalContractorLocation(assignee_key);

        if( assignee_location == Location::NoWhere() )
        {
            std::cout<<__func__<<": cannot find assignee! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        mAssignorInputs[assignee_location].push_back(assignment_input);

        return assignment_key;
    }

    void ExecuteAllDistributedAssignments()
    {
        SendAssignorInputsToAssigneeInputs(0);
        GenerateWorkUnitsInputs();
        ExecuteAllLocalWorkUnits();
        GenerateAssigneeOuputs();
        SendAssigneeOutputsToAssignorOutputs(1);
    }

    void SendAssignorInputsToAssigneeInputs( const int mpi_tag = 0 )
    {
        mpCommunicator->AllSendAllRecv( mAssignorInputs, mAssigneeInputs, mpi_tag );
    }

    void GenerateWorkUnitsInputs()//convert assignee input data to work unit input data
    {
        typedef AssignmentDataType<TInputType> Input;
        typedef AssignmentDataVectorType<TInputType> InputVector;
        typedef AssignmentDataVectorPairType<Location, TInputType> InputVectorPairByLocation;

        mWorkUnitInputs.clear();

        for( const InputVectorPairByLocation & r_assignee_input_vector_pair : mAssigneeInputs )
        {
            const InputVector & r_assignee_input_vector = r_assignee_input_vector_pair.second;
            for ( const Input & r_assignee_input : r_assignee_input_vector )
            {
                ContractorKey assignee_key = r_assignee_input.mAssigneeKey;
                mWorkUnitInputs[assignee_key].push_back(r_assignee_input);
            }
        }
    }

    void ExecuteAllLocalWorkUnits()
    {
        using InputVectorPairByContractorKey = AssignmentDataVectorPairType<ContractorKey, TInputType> ;

        for( const InputVectorPairByContractorKey & r_input_vector_pair : mWorkUnitInputs )
        {
            const ContractorKey assignee_key = r_input_vector_pair.first;
            ExecuteLocalWorkUnit(assignee_key);
        }
    }

    void ExecuteLocalWorkUnit( const ContractorKey assignee_key )
    {
        typedef AssignmentDataVectorType<TInputType> InputVector;
        typedef AssignmentDataVectorType<TOutputType> OutputVector;
        typedef AssignmentDataVectorMapType<ContractorKey, TInputType>  InputVectorMap;
        typedef typename InputVectorMap::const_iterator InputVectorMapIterator;

        //assignee
        TAssigneeType * p_assignee = mpAssigneeManager->FindLocalContractorPointer(assignee_key);
        if( p_assignee == nullptr )
        {
            std::cout <<__func__<<"local contractor not found! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }

        //work unit inputs
        InputVectorMapIterator it_input_vector_map = mWorkUnitInputs.find(assignee_key);
        if( it_input_vector_map == mWorkUnitInputs.end() )
        {
            std::cout<<__func__<<": work unit inputs not found! exit"<<std::endl;
            exit(EXIT_FAILURE);
        }
        const InputVector & r_input_vector = it_input_vector_map->second;
        typename  InputVector::size_type input_vector_size = r_input_vector.size();

        //work unit outputs
        OutputVector & r_output_vector = mWorkUnitOutputs[assignee_key];
        typename OutputVector::size_type output_vector_size = (typename OutputVector::size_type) input_vector_size;

        //clear work unit outputs
        r_output_vector.clear();
        r_output_vector.resize(output_vector_size);

        //execute
        {
            typename InputVector ::size_type i_input;
            typename OutputVector::size_type i_output;

            for ( i_input = 0, i_output = 0; i_input < input_vector_size || i_output < output_vector_size; i_input++, i_output++ )
            {
                p_assignee->Execute(r_input_vector[i_input].mData, r_output_vector[i_output].mData);

                r_output_vector[i_output].mAssignorKey   = r_input_vector[i_input].mAssignorKey;
                r_output_vector[i_output].mAssignmentKey = r_input_vector[i_input].mAssignmentKey;
                r_output_vector[i_output].mAssigneeKey   = r_input_vector[i_input].mAssigneeKey;
            }
        }
    }

    void GenerateAssigneeOuputs()//convert work unit output data to assignee output data
    {
        typedef AssignmentDataType<TOutputType> Output;
        typedef AssignmentDataVectorType<TOutputType> OutputVector;
        typedef AssignmentDataVectorPairType<ContractorKey, TOutputType> OutputVectorPairByContractorKey;

        mAssigneeOutputs.clear();

        for( const OutputVectorPairByContractorKey & r_work_unit_output_vector_pair : mWorkUnitOutputs )
        {
            const OutputVector & r_work_unit_output_vector = r_work_unit_output_vector_pair.second;
            for ( const Output & r_output : r_work_unit_output_vector )
            {
                ContractorKey assignor_key = r_output.mAssignorKey;
                Location assignor_location = mpAssignorManager->FindGlobalContractorLocation(assignor_key);

                if( assignor_location == Location::NoWhere() )
                {
                    std::cout<<__func__<<": cannot find assignor! exit"<<std::endl;
                    exit(EXIT_FAILURE);
                }

                mAssigneeOutputs[assignor_location].push_back(r_output);
            }
        }
    }

    void SendAssigneeOutputsToAssignorOutputs( const int mpi_tag = 0 )
    {
        mpCommunicator->AllSendAllRecv( mAssigneeOutputs, mAssignorOutputs, mpi_tag );
    }

    void GetResultsAtAssignor( AssignmentDataVectorType<TOutputType> & r_result_vector ) const
    {
        using Output = AssignmentDataType<TOutputType> ;
        using OutputVector = AssignmentDataVectorType<TOutputType> ;
        using OutputVectorPairByLocation = AssignmentDataVectorPairType<Location, TOutputType> ;

        r_result_vector.clear();

        for( const OutputVectorPairByLocation & r_output_vector_pair : mAssignorOutputs )
        {
            const OutputVector & r_output_vector = r_output_vector_pair.second;
            for( const Output & r_output : r_output_vector )
                r_result_vector.push_back(r_output);
        }
    }

    void GetResultsAtAssignee( AssignmentDataVectorType<TOutputType> & r_result_vector ) const
    {
        using Output = AssignmentDataType<TOutputType> ;
        using OutputVector = AssignmentDataVectorType<TOutputType> ;
        using OutputVectorPairByContractorKey = AssignmentDataVectorPairType<ContractorKey, TOutputType> ;

        r_result_vector.clear();

        for( const OutputVectorPairByContractorKey & r_output_vector_pair : mWorkUnitOutputs )
        {
            const OutputVector & r_output_vector = r_output_vector_pair.second;
            for( const Output & r_output : r_output_vector )
                r_result_vector.push_back(r_output);
        }
    }

    void PrintAllAssignments() const
    {
        std::cout << __func__ << ": Assignor input vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<Location, TInputType> & r_assignor_input_vector_pair : mAssignorInputs )
        {
            const Location & r_location = r_assignor_input_vector_pair.first;
            const AssignmentDataVectorType<TInputType> & r_assignor_input_vector = r_assignor_input_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Location: ";
            printer.Print(r_location);
            std::cout<<"Assignor input vector: ";
            printer.Print(r_assignor_input_vector);
        }

        std::cout << std::endl;
        std::cout << __func__ << ": Assignee input vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<Location, TInputType> & r_assignee_input_vector_pair : mAssigneeInputs )
        {
            const Location & r_location = r_assignee_input_vector_pair.first;
            const AssignmentDataVectorType<TInputType> & r_assignee_input_vector = r_assignee_input_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Location: ";
            printer.Print(r_location);
            std::cout<<"Assignee input vector: ";
            printer.Print(r_assignee_input_vector);
        }

        std::cout << std::endl;
        std::cout << __func__ << ": Work Unit input vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<ContractorKey, TInputType> & r_work_unit_input_vector_pair : mWorkUnitInputs )
        {
            const ContractorKey & r_assignee_key = r_work_unit_input_vector_pair.first;
            const AssignmentDataVectorType<TInputType> & r_work_unit_input_vector = r_work_unit_input_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Assignee Key: ";
            printer.Print(r_assignee_key);
            std::cout<<"Assignee input vector: ";
            printer.Print(r_work_unit_input_vector);
        }

        std::cout << std::endl;
        std::cout << __func__ << ": Work Unit output vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<ContractorKey, TOutputType> & r_work_unit_output_vector_pair : mWorkUnitOutputs )
        {
            const ContractorKey & r_assignee_key = r_work_unit_output_vector_pair.first;
            const AssignmentDataVectorType<TOutputType> & r_work_unit_output_vector = r_work_unit_output_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Assignee Key: ";
            printer.Print(r_assignee_key);
            std::cout<<"Assignee output vector: ";
            printer.Print(r_work_unit_output_vector);
        }

        std::cout << std::endl;
        std::cout << __func__ << ": Assignee output vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<Location, TOutputType> & r_assignee_output_vector_pair : mAssigneeOutputs )
        {
            const Location & r_location = r_assignee_output_vector_pair.first;
            const AssignmentDataVectorType<TOutputType> & r_assignee_output_vector = r_assignee_output_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Location: ";
            printer.Print(r_location);
            std::cout<<"Assignee output vector: ";
            printer.Print(r_assignee_output_vector);
        }

        std::cout << std::endl;
        std::cout << __func__ << ": Assignor output vector map: " << std::endl;

        for( const AssignmentDataVectorPairType<Location, TOutputType> & r_assignor_output_vector_pair : mAssignorOutputs )
        {
            const Location & r_location = r_assignor_output_vector_pair.first;
            const AssignmentDataVectorType<TOutputType> & r_assignor_output_vector = r_assignor_output_vector_pair.second;

            DataUtility::DataPrinter printer;
            std::cout<<"Location: ";
            printer.Print(r_location);
            std::cout<<"Assignor output vector: ";
            printer.Print(r_assignor_output_vector);
        }

        std::cout << std::endl;
    }
private:
    TCommunicatorType * mpCommunicator;

    ContractorManagerType<TAssignorType> * mpAssignorManager;
    ContractorManagerType<TAssigneeType> * mpAssigneeManager;

    AssignmentKeyIssuer mAssignmentKeyIssuer;

    AssignmentDataVectorMapType<Location,           TInputType> mAssignorInputs;
    AssignmentDataVectorMapType<Location,           TInputType> mAssigneeInputs;
    AssignmentDataVectorMapType<ContractorKey, TInputType> mWorkUnitInputs;

    AssignmentDataVectorMapType<Location,           TOutputType> mAssignorOutputs;
    AssignmentDataVectorMapType<Location,           TOutputType> mAssigneeOutputs;
    AssignmentDataVectorMapType<ContractorKey, TOutputType> mWorkUnitOutputs;
};

}
